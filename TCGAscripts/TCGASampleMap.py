import string, os, sys

sys.path.insert(0,"../CGDataNew")

from SampleMapNew import *
from CGDataLib import *
from ClinicalMatrixNew import *
from IntegrationId import *
import TCGAUtil

def TCGASampleMap (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, __name__

    #if cancer in ["PANCAN","PANCAN12"]:
    #    return

    ignore =1
    bookDic = cgWalk(dir,ignore)
    
    existMaps = collectSampleMaps(bookDic)
    missingMaps=  collectMissingSampleMaps(bookDic)

    #removeExistMaps
    for map in existMaps:
        if map not in missingMaps:
            missingMaps[map]=existMaps[map]
        
    # all aliquote uuid dic
    aliquote_dic =TCGAUtil.uuid_Aliquot_all()
    sample_dic =TCGAUtil.uuid_Sample_all()
    
    #missingMaps --- actually this is all the maps
    for map in missingMaps:
        print map
        print missingMaps[map]
        sMap =SampleMapNew(None,map)

        #integration id
        intName= map+".integrationID"
        if intName in bookDic:
            fin = open(bookDic[intName]["path"],'r')
            integrationID=IntegrationId(intName, fin)
            fin.close()
        else:
            integrationID=IntegrationId(intName)

        samples =[]
        for name in missingMaps[map]:
            if REALRUN !=1:
                continue
            print name
            obj=bookDic[name]
            if obj['type']=="genomicMatrix":
                fin =open(obj['path'],'U')
                for sample in string.split(fin.readline()[:-1],"\t")[1:]:
                    if sample =="":
                        print name, "has bad empty sample id"
                    if sample not in samples:
                        samples.append(sample)
                fin.close()
            #elif obj['type']=="clinicalMatrix":
            #    cMa = ClinicalMatrixNew(obj['path'],name)
            #    for sample in cMa.getROWs():
            #        if sample not in samples:
            #            samples.append(sample)
            elif obj['type'] in ["mutationVector","clinicalMatrix"]:
                path = obj['path']
                os.system("cut -f 1 "+path+ " |sort |uniq > .tmp")
                fin=open('.tmp','r')
                fin.readline()
                for line in fin.readlines():
                    #if string.strip(line)=="":
                    #    break
                    sample = string.strip(line) #string.split(line,'\t')[0]
                    if sample =="":
                        break
                    if sample not in samples:
                        samples.append(sample)

            else:
                continue

        for sample in samples:
            if REALRUN !=1:
                continue
            #TCGA uuid handling
            if sample[0:4]!="TCGA": 
                if aliquote_dic.has_key(string.lower(sample)):
                    TCGAbarcode = aliquote_dic[string.lower(sample)]
                else:
                    print sample
                    continue
                parent = TCGAbarcode
                child = sample
                sMap.addLink(parent,string.lower(child))
                sMap.addLink(parent,string.upper(child))
                sample = parent

            #do TCGA barcode trick
            parts= string.split(sample,"-")
            if len(parts)>3 and len(parts[3])==3:
                parts = parts[0:3]+ [parts[3][0:2],parts[3][2]]+parts[4:]
                #print parts

            """
            parent = string.join(parts[0:3],"-")
            #parts[3]
            if len(parts)>3 and len(parts[3])==3:
                child=parent +"-" +parts[3][0:2]
                sMap.addLink(parent,child)
                parent=child
                child=string.join(parts[0:4],"-")
                sMap.addLink(parent,child)
                parent=child
            """
            parent = string.join(parts[0:3],"-")
            for i in range (3,len(parts)):
                if i!=4:
                    child = parent +"-" +parts[i]
                else:
                    child = parent +parts[i]
                #add parent child
                sMap.addLink(parent,child)
                parent = child
                
            intID= TCGAUtil.barcode_IntegrationId(sample)
            integrationID.addId(intID)
            
        #output sampleMap
        if not os.path.exists( outDir ):
            os.makedirs( outDir )
        if not os.path.exists( outDir +cancer+"/"):
                os.makedirs( outDir+cancer+"/" )

        if REALRUN == 1:
            oHandle = open(outDir+cancer+"/"+map,"w")
            sMap.store(oHandle)

        #output integrationID
        if REALRUN ==1:
            oHandle = open(outDir+cancer+"/integrationID","w")
            integrationID.store(oHandle)
            oHandle.close()
        
        #output integrationID json
        oHandle = open(outDir+cancer+"/integrationID.json","w")
        J={}
        J['name']=intName

        J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
        J["sample_type"]="tumor"
        if cancer not in ["PANCAN","PANCAN12"]:
            J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
        else:
            J["primary_disease"]="cancer"
            
        #J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]
        J['domain']="TCGA"
        J['owner']="TCGA"
        
        J["cgDataVersion"]=1
        J['type']="integrationId"
        J["version"]= datetime.date.today().isoformat()
        oHandle.write( json.dumps( J, indent=-1 ) )
        oHandle.close()
        
        #output json
        oHandle = open(outDir+cancer+"/"+map+".json","w")
        J['name']=map
        J['type']="sampleMap"
        J["version"]= datetime.date.today().isoformat()
        J["cgDataVersion"]=1
        J[":integrationId"]=intName

        #add info for old clinical data
        if os.path.exists( outDir+cancer+"/oldClin.json" ):
            J[':oldClin']=cancer+"_oldClin" 

        #special code
        if TCGAUtil.featurePriority.has_key(cancer) and len(TCGAUtil.featurePriority[cancer])>=5:
            J["VIS"]=5
        
        #blackList in PAAD
        if J['name'] in ["TCGA.PAAD.sampleMap"]:
            J["blacklist"]= [ "TCGA-FQ-6551",
                              "TCGA-FQ-6552",
                              "TCGA-FQ-6553",
                              "TCGA-FQ-6554",
                              "TCGA-FQ-6555",
                              "TCGA-FQ-6558",
                              "TCGA-FQ-6559"]
            
        oHandle.write( json.dumps( J, indent=-1 ) )

        
    return
