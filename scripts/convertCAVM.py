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
        if string.find(os.path.abspath(path), "/inside/home/jzhu/cgDataJing/scripts/data_flatten/") ==-1:
            print "ignore "+path
            continue
            
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

        for name in sampleMaps[sampleMap]:
            obj=bookDic[name]
            if obj['type']=="clinicalMatrix":
                if REALRUN ==-1:
                    continue
                
                outfile = outDir +os.path.basename(obj['path'])
                fin = open(obj['path'],'r')
                fout = open(outfile,'w')
                fout.write(fin.readline())

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
                integrationList = samples[:]
                
                #JSON
                fin = open (obj['path']+".json",'r')
                J=json.load(fin)
                fin.close()
                J['cohort']=J[':sampleMap']
                if CAVM:
                    J.pop(':sampleMap') 
                fout=open(outfile+".json",'w')
                fout.write(json.dumps (J, indent=-1))
                fout.close()

                #clinicalFeature
                if J.has_key(":clinicalFeature"):
                    cFobj= bookDic[J[":clinicalFeature"]]
                    outfile = outDir +os.path.basename(cFobj['path'])
                    os.system("cp "+cFobj['path']+"  "+outfile)
                    os.system("cp "+cFobj['path']+".json "+outfile+".json")

                
                #sampleMap data #cgData 1
                if not CAVM:
                    outfile = outDir +"sampleMap"
                    fout=open(outfile,'w')
                    for sample in samples:
                        fout.write(sample+"\t"+sample+"\n")
                    fout.close()
                
        for name in sampleMaps[sampleMap]:
            obj=bookDic[name]
            if obj['type']=="genomicSegment":
                path= obj['path']
                os.system("cp "+path+" "+outDir+os.path.basename(path))
                os.system("cp "+path+".json "+outDir+os.path.basename(path)+".json")
                
            if obj['type']=="genomicMatrix":
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
                    root = sMap.getIntegrationId(sample,integrationList)
                    if roots.has_key(root):
                        roots[root].append(i)
                        findDup=1
                    else:
                        roots[root]=[i]
                fin.close()

                #JSON
                outfile = outDir +os.path.basename(obj['path'])
                fin = open (obj['path']+".json",'r')
                J=json.load(fin)
                fin.close()
                J['cohort']=J[':sampleMap']
                if CAVM:
                    J.pop(':sampleMap')
                fout=open(outfile+".json",'w')
                fout.write(json.dumps (J, indent=-1))
                fout.close()

                if J.has_key('anatomical_origin'):
                    sMapJ['anatomical_origin']=J['anatomical_origin']
                    if CAVM:
                        J.pop('anatomical_origin')
                if J.has_key('primary_disease'):
                    sMapJ['primary_disease']=J['primary_disease']
                    if CAVM:
                        J.pop('primary_disease') 
                if J.has_key('domain'):
                    sMapJ['domain']=J['domain']
                    if CAVM:
                        J.pop('domain') 
                if J.has_key('sample_type'):
                    sMapJ['sample_type']=J['sample_type']
                    if CAVM:
                        J.pop('sample_type')
                if J.has_key('tags'):
                    sMapJ['tags']=J['tags']
                    if CAVM:
                        J.pop('tags')

                if J.has_key(':genomicSegment'):
                    probeMap = bookDic[J[':probeMap']]['path']
                    os.system("cp "+probeMap+" "+outDir+os.path.basename(probeMap))
                    os.system("cp "+probeMap+".json "+outDir+os.path.basename(probeMap)+".json")
                    
                if REALRUN != 1:
                    continue

                fin =open(obj['path'],'U')
                fout=open(outfile, 'w')
                #genomic data no dup
                if findDup ==0:
                    data =string.split(fin.readline()[:-1],"\t")
                    samples= data[1:]
                    fout.write(data[0])
                    for sample in samples:
                        fout.write('\t'+sMap.getIntegrationId(sample,integrationList))
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
                    while 1:
                        line = fin.readline()[:-1]
                        if line =="":
                            break
                        data = string.split(line,"\t")
                        fout.write(data[0])
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
                            fout.write('\t'+average)
                        fout.write('\n')
                    fout.close()
                    
                
        #sampleMap json #cgData1
        if not CAVM:
            outfile = outDir +"sampleMap.json"
            fout=open(outfile, 'w')
            fout.write(json.dumps(sMapJ,indent=-1))
            fout.close()
        
        #cohort json
        outfile= outDir+"cohort.json"
        fout=open(outfile,'w')
        cohortJ = copy.deepcopy(sMapJ)
        cohortJ['type']="cohort"
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
    print "          --REALRUN=0 (full)"
    print "          --REALRUN=1 (clinical data probeMap genomicSegment only, full json, skip genomic data )"
    print "          --REALRUN=-1 (probeMap genomicSegment only, full genomic json, skip genomic data )"
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


