import string, os,sys, json
import optparse

sys.path.insert(0,"../CGDataNew")

from SampleMapNew import *
from CGDataLib import *

def convertCAVM (inDir, outD ,REALRUN, CAVM, TCGA, MAPID=1):
    bookDic={}
    sampleMaps={}
    bookDic=cgWalk(inDir,0)
    
    if not os.path.exists (outD):
        os.system("mkdir "+outD)

    if not bookDic :
        print "repo has problem"
        return 0


    sampleMaps = collectSampleMaps(bookDic)
    missingMaps= collectMissingSampleMaps(bookDic)

    allMaps = sampleMaps.keys()
    allMaps.extend(missingMaps.keys())

    for sampleMap in allMaps:
        print sampleMap
        outDir = outD + sampleMap+"/"
        if not os.path.exists (outDir):
            os.system("mkdir "+outDir)

        path = bookDic[sampleMap]['path']
        """
        if string.find(os.path.abspath(path), "/inside/home/jzhu/cgDataJing/scripts/data_flatten/") ==-1:
            print "ignore "+path
            continue
        """    
        if sampleMap in missingMaps:
            #construct an empty sampleMap
            sMap = SampleMapNew(None,sampleMap)
            #fill sMap with individual nodes, no connection
            changed = checkIdsAllIn(sMap, bookDic)
            #build connection
        else:
            name = bookDic[sampleMap]['name']
            fin = open(path,'r')
            sMap = SampleMapNew(fin,name)
            if not sMap.getName():
                print "Fail to initiate", name
                return 0
            fin.close()

        #cohort sampleMap json
        sMapJ={}
        fin = open(bookDic[sampleMap]['path']+".json",'r')
        sMapJ = json.loads(fin.read())
        fin.close()

        #integration list
        integrationList =[]
        rootDic={}
        clinFile =""
        clinMatrix = None
        
        #cohort
        COHORT =""
        cohortPath= string.join(string.split(bookDic[sampleMap]['path'],"/")[0:-1],"/")+"/cohort.json"
        if os.path.exists(cohortPath):
            fin = open (cohortPath,'r')
            cohortJ= json.loads(fin.read())
            COHORT = cohortJ["name"]

        for name in sampleMaps[sampleMap]:
            obj=bookDic[name]
            if obj['type']=="clinicalMatrix":
                clinFile = outDir +os.path.basename(obj['path'])

                #JSON
                fin = open (obj['path']+".json",'r')
                J=json.load(fin)
                fin.close()

                if COHORT :
                    J["cohort"]=COHORT
                else:
                    J['cohort']=J[':sampleMap']

                J["label"]="Phenotypes"
                if CAVM:
                    J.pop(':sampleMap') 
                    if J.has_key("dataSubType"):
                        if J.has_key(":dataSubType"):
                            J.pop(':dataSubType') 
                    else:
                        if J.has_key(":dataSubType"):
                            J["dataSubType"]=J[":dataSubType"]
                            J.pop(':dataSubType')

                fout=open(clinFile+".json",'w')
                fout.write(json.dumps (J, indent=-1))
                fout.close()

                if REALRUN != 0 and  REALRUN !=1:
                    continue

                if clinMatrix != None:
                    print "only one clinical matrix is allowed"
                    sys.exit()
                    
                fin = open(obj['path'],'U')
                fout = open(clinFile,'w')
                line = fin.readline()
                fout.write(line)

                samples =[]
                for line in fin.readlines():
                    if MAPID:
                        sample =string.split(line[:-1],"\t")[-1]
                    else:
                        sample =string.split(line[:-1],"\t")[0]
                    if sample not in samples and sample !="":
                        samples.append(sample)
                        fout.write(sample+"\t")
                        if MAPID:
                            fout.write(string.join(string.split(line[:-1],"\t")[1:],"\t"))
                        else:
                            fout.write(string.join(string.split(line[:-1],"\t")[1:-1],"\t")+"\t"+sample)
                        fout.write("\n")
                fout.close()
                integrationList = copy.deepcopy(samples)
                

                #clinicalFeature
                if J.has_key(":clinicalFeature"):
                    cFobj= bookDic[J[":clinicalFeature"]]
                    outfile = outDir +os.path.basename(cFobj['path'])
                    os.system("cp "+cFobj['path']+"  "+outfile)
                    os.system("cp "+cFobj['path']+".json "+outfile+".json")
                
                #sampleMap data mapping information #cgData 1
                if not CAVM:
                    os.system("cp "+ bookDic[sampleMap]['path'] +" "+outDir+"sampleMap")

                #only expect one clinical matrix
                clinMatrix= ClinicalMatrixNew(clinFile, "clinMatrix")
                break

        for name in sampleMaps[sampleMap]:
            obj=bookDic[name]
            if obj['type'] in ["genomicSegment","mutationVector"]:
                path= obj['path']
                print path

                outfile = outDir +os.path.basename(obj['path'])
                fin = open (obj['path']+".json",'r')
                J=json.load(fin)
                fin.close()

                if COHORT :
                    J["cohort"]=COHORT
                else:
                    J['cohort']=J[':sampleMap']

                if CAVM:
                    J.pop(':sampleMap')
                    if J.has_key("dataSubType"):
                        if J.has_key(":dataSubType"):
                            J.pop(':dataSubType') 
                    else:
                        if J.has_key(":dataSubType"):
                            J["dataSubType"]=J[":dataSubType"]
                            J.pop(':dataSubType') 


                fout=open(outfile+".json",'w')
                fout.write(json.dumps (J, indent=-1))
                fout.close()

                if REALRUN ==1 :
                    fin =open(path,'r')
                    fout =open(outDir+os.path.basename(path),'w')
                    for line in fin.readlines():
                        data = string.split(line,'\t')
                        sample =data[0]
                        if rootDic.has_key(sample):
                            root = rootDic[sample]
                        else:
                            root = sMap.getIntegrationId(sample,integrationList)
                            if not root:
                                root = sample
                            rootDic[sample]=root
                        fout.write(root+"\t"+string.join(data[1:],'\t'))
                    fin.close()
                    fout.close()
                
            if obj['type']=="genomicMatrix":
                print obj['name']
                #JSON
                outfile = outDir +os.path.basename(obj['path'])
                fin = open (obj['path']+".json",'r')
                J=json.load(fin)
                fin.close()

                if COHORT :
                    J["cohort"]=COHORT
                else:
                    J['cohort']=J[':sampleMap']

                if CAVM:
                    J.pop(':sampleMap')
                    if J.has_key("dataSubType"):
                        if J.has_key(":dataSubType"):
                            J.pop(':dataSubType') 
                    else:
                        if J.has_key(":dataSubType"):
                            J["dataSubType"]=J[":dataSubType"]
                            J.pop(':dataSubType') 

                fout=open(outfile+".json",'w')
                fout.write(json.dumps (J, indent=-1))
                fout.close()

                if J.has_key('anatomical_origin'):
                    sMapJ['anatomical_origin']=J['anatomical_origin']
                if J.has_key('primary_disease'):
                    sMapJ['primary_disease']=J['primary_disease']
                if J.has_key('domain'):
                    sMapJ['domain']=J['domain']
                if J.has_key('sample_type'):
                    sMapJ['sample_type']=J['sample_type']
                if J.has_key('tags'):
                    sMapJ['tags']=J['tags']
                
                if REALRUN != 1 and REALRUN !=0:
                    continue
                
                # add to clinMatrix the id mappings
                mappingCol= "_GENOMIC_ID_"+obj['name']
                clinMatrix.addOneColWithSameValue(mappingCol,"")
                
                # need to find it out if there are more than one sample map to each _INTEGRATION ID
                roots={}
                findDup=0
                fin =open(obj['path'],'U')
                samples =string.split(fin.readline()[:-1],"\t")[1:]
                for i in range(0,len(samples)):
                    sample = samples[i]
                    if sample =="":
                        print name, "has bad empty sample id"
                        sys.exit()
                    if rootDic.has_key(sample):
                        root = rootDic[sample]
                    else:
                        root = sMap.getIntegrationId(sample,integrationList)
                        if not root:
                            root = sample
                        rootDic[sample]=root

                    genomic_Id= clinMatrix.getDATA(root, mappingCol)
                    if genomic_Id is None or genomic_Id =="":
                        clinMatrix.setDATA(root, mappingCol,sample)
                    else:
                        genomic_Id = string.split(genomic_Id,",")
                        if sample not in genomic_Id:
                            genomic_Id.append(sample)
                            genomic_Id= string.join(genomic_Id,',')
                            #print sample, genomic_Id
                            clinMatrix.setDATA(root, mappingCol,genomic_Id)
                            
                    if roots.has_key(root):
                        roots[root].append(i)
                        findDup=1
                    else:
                        roots[root]=[i]
                fin.close()

                if REALRUN != 1:
                    continue
                
                #probemap for genomic segment
                if J.has_key(':genomicSegment'):
                    probeMap = bookDic[J[':probeMap']]['path']
                    os.system("cp "+probeMap+" "+outDir+os.path.basename(probeMap))
                    os.system("cp "+probeMap+".json "+outDir+os.path.basename(probeMap)+".json")

                #need to figure out if there are duplication in the probe ids
                findDupProbe=[]
                process = os.popen("r=$(cut -f 1  "+obj['path']+" | more +2 | sort |uniq -c | sed -e 's/ *//' -e 's/ /\t/' | sort -n |cut -f 1 |sort -un|tail -n 1); if [ $r -ne \"1\" ]; then echo $r ; fi")
                r = process.read()
                if r:
                    print string.strip(r), obj['path']
                    process = os.popen("cut -f 1  "+obj['path']+" | more +2 | sort |uniq -c | sed -e 's/ *//' -e 's/ /\t/' | sort -n |  grep -vP ^1'\t' | cut -f 2 |sort")
                    r = process.read()
                    list = string.split(r,"\n")
                    print len(list)
                    for probe in list:
                        findDupProbe.append(probe)
                        
                #genomic data no dup
                fout=open(outfile ,'w')
                fin =open(obj['path'],'U')
                if findDup ==0 and findDupProbe ==[]:
                    data =string.split(fin.readline()[:-1],"\t")
                    samples= data[1:]
                    fout.write(data[0])
                    for sample in samples:
                        if rootDic.has_key(sample):
                            root = rootDic[sample]
                        else:
                            root = sMap.getIntegrationId(sample,integrationList)
                            if not root:
                                root = sample
                            rootDic[sample]=root
                        fout.write('\t'+root)
                    fout.write('\n')

                    if TCGA:
                        fin.close()
                        fout.close()
                        os.system("more +2 "+obj['path']+" >> "+outfile)
                    else:
                        while 1:
                            line = fin.readline()
                            if line =="":
                                break
                            line = string.replace(line,"\tnan\t","\tNA\t")
                            line = string.replace(line,"\tNAN\t","\tNA\t")
                            line = string.replace(line,"\tNaN\t","\tNA\t")
                            fout.write(line)
                        fin.close()
                        fout.close()

                #genomic data with dup
                else:
                    print "genomic with dup",obj['path']
                    data =string.split(fin.readline()[:-1],"\t")
                    fout.write(data[0])
                    for root in roots:
                        fout.write('\t'+root)
                    fout.write('\n')

                    dupDic ={}
                    while 1:
                        duplist=[]
                        line = fin.readline()[:-1]
                        if line =="":
                            break

                        data = string.split(line,"\t")
                                                
                        if data[0] not in findDupProbe:
                            fout.write(data[0])                            
                        else:
                            if data[0] not in dupDic:
                                dupDic[data[0]]=[]

                        values =data[1:]
                        
                        for root in roots:
                            if len(roots[root])!=1:
                                total="NA"
                                n=0
                                for i in roots[root]:
                                    if values[i] in ["nan","NAN","NaN"]:
                                        pass
                                    else:
                                        try:
                                            float(values[i])
                                            if total=="NA":
                                                total = float(values[i])
                                            else:
                                                total = total +float(values[i])
                                            n=n+1
                                        except:
                                            pass
                                if total != "NA":
                                    average = str(total / n)
                                else:
                                    average ="NA"
                            else:
                                if values[roots[root][0]] in ["nan","NAN","NaN"]:
                                    average="NA"
                                else:
                                    try:
                                        float(values[roots[root][0]])
                                        average = values[roots[root][0]]
                                    except:
                                        average="NA"
                            if data[0] not in findDupProbe:
                                fout.write('\t'+average)
                            else:
                                duplist.append(average)
                        if data[0] not in findDupProbe:
                            fout.write('\n')
                        else:
                            dupDic[data[0]].append(duplist[:])
                            
                    if dupDic!={}:
                        for probe in dupDic:
                            fout.write(probe)
                            k = len (dupDic[probe][0])
                            valList=[]
                            nList=[]
                            for i in range (0,k):
                                valList.append("NA")
                                nList.append(0)
                            
                            for list in dupDic[probe]:
                                for i in range (0,k):
                                    try:
                                        float(list[i])
                                        if valList [i]=="NA":
                                            valList[i]=float(list[i])
                                        else:
                                            valList[i] =valList[i] +float(list[i])
                                        nList[i]=nList[i]+1
                                    except:
                                        pass
                            for i in range (0,k):
                                try:
                                    float(valList[i])
                                    fout.write("\t"+str(float(valList[i])/nList[i]))
                                except:
                                    fout.write("\tNA")
                            fout.write("\n")
                    fin.close()
                    fout.close()
                

        #final clinical matrix output
        if REALRUN == 0 or  REALRUN ==1:
            fout= open(clinFile,'w')
            clinMatrix.store(fout)
    
        #sampleMap json #cgData1
        if not CAVM:
            outfile = outDir +"sampleMap.json"
            fout=open(outfile, 'w')
            fout.write(json.dumps(sMapJ,indent=-1))
            fout.close()
        
        #cohort json cp or create
        cohortPath= string.join(string.split(bookDic[sampleMap]['path'],"/")[0:-1],"/")+"/cohort.json"
        if os.path.exists(cohortPath):
            os.system("cp " + cohortPath +" " + outDir)
        else:
            outfile = outDir+"cohort.json"
            fout=open(outfile,'w')
            cohortJ={}
            cohortJ["type"]="cohort"
            cohortJ["name"]= sampleMap
            fout.write(json.dumps(cohortJ,indent=-1))
            fout.close()

parser = optparse.OptionParser()
parser.add_option("--inDir", action="store", type="string", dest="inDir")
parser.add_option("--outDir", action="store", type="string", dest="outDir")
parser.add_option("--CAVM", action="store", type="string", dest="CAVM")
parser.add_option("--TCGA", action="store", type="string", dest="TCGA")
parser.add_option("--REALRUN", action="store", type="string", dest="REALRUN")
parser.add_option("--MAPID", action="store", type="string", dest="MAPID")
(options, args) = parser.parse_args()

def printUsage():
    print "python convertCAVM.py --inDir=inputDir --outDir=outputDir --CAVM=0/1 --TCGA=0/1 --REALRUN=0/1/-1"
    print
    print "          --CAVM=0 (for old stype cgData)"
    print "          --CAVM=1 (for CAVM stype cgData)"
    print "          --TCGA=0 (input data is not TCGA data)"
    print "          --TCGA=1 (intput data is TCGA data, skip checking NaNs)"
    print "          --REALRUN=1 (full)"
    print "          --REALRUN=0 (clinical data probeMap only, full json, skip genomic data, genomic segment data)"
    print "          --REALRUN=-1 (full json, skip all genomic data and clinical data)"
    print
    print " options: "
    print "          --MAPID=0 (do not change genomic id)"
    print 


if options.inDir==None or options.outDir ==None \
       or options.CAVM ==None or options.TCGA==None \
       or options.REALRUN ==None:
    printUsage()
    sys.exit()


inDir = options.inDir
outDir =options.outDir

if inDir[-1]!="/":
    inDir = inDir +"/"

if outDir[-1]!="/":
    outDir = outDir +"/"

REALRUN=int(options.REALRUN)
CAVM =int (options.CAVM)
TCGA = int(options.TCGA)
if options.MAPID==None:
    MAPID= 1
else:
    MAPID= int(options.MAPID)

convertCAVM (inDir,outDir, REALRUN, CAVM, TCGA, MAPID)


