import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def cohort (inDir, outDir, cancer, flog,REALRUN):
    #print status
    if cancer in ["COADREAD","LUNG","PANCAN"]:
        return

    print cancer, sys._getframe().f_code.co_name
        
    print inDir
    print outDir

    if REALRUN:
        ignore =1
        bookDic=cgWalk(inDir,ignore)
        
        existMaps = collectSampleMaps(bookDic)
        missingMaps=  collectMissingSampleMaps(bookDic)

        #removeExistMaps
        for map in existMaps:
            if map not in missingMaps:
                missingMaps[map]=existMaps[map]
        
        # all aliquote uuid dic
        aliquote_dic =TCGAUtil.uuid_Aliquot_all()
        sample_dic =TCGAUtil.uuid_Sample_all()

        if len(missingMaps)!=1:
            return

        map = missingMaps.keys()[0]
        print map
        samples =[]
        for name in missingMaps[map]:
            obj=bookDic[name]
                
            if obj['type']=="genomicMatrix":
                fin =open(obj['path'],'U')
                for sample in string.split(fin.readline()[:-1],"\t")[1:]:
                    if sample =="":
                        print name, "has bad empty sample id"
                        sys.exit()
                    if sample not in samples:
                        samples.append(sample)
                fin.close()
                
        intDic={}
        for sample in samples:
            #TCGA uuid handling
            uuid =sample
            TCGAbarcode =""
            if uuid[0:4]!="TCGA": 
                if aliquote_dic.has_key(string.lower(uuid)):
                    TCGAbarcode = aliquote_dic[string.lower(uuid)]
                else:
                    TCGAbarcode =  uuid
            else:
                TCGAbarcode = sample

            intID= TCGAUtil.barcode_IntegrationId(TCGAbarcode)
            if intID == None: # ids is on patient level above integration level
                continue 
            if not intDic.has_key(intID):
                intDic[intID]=""

        outfile = outDir+cancer+"/cohort"
        fout =open(outfile,"w")
        fout.write("sample\tcohort\n")
        for intId in intDic:
            fout.write(intId+"\t"+TCGAUtil.cancerHumanReadable[cancer]+"\n")
        fout.close()

    #clinFeature
    feature="cohort"
    withClinF =0
    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(feature):
            priority = TCGAUtil.featurePriority[cancer][feature]

            J={}
            J["name"]="TCGA_"+cancer+"_cohort_clinFeat"
            J["type"]="clinicalFeature"
    
            oHandle = open(outDir+cancer+"/cohort_clinFeature.json","w")
            oHandle.write( json.dumps( J, indent=-1 ) )
            oHandle.close()

            oHandle = open(outDir+cancer+"/cohort_clinFeature","w")
            oHandle.write(feature+"\tpriority\t"+str(priority)+"\n")
            oHandle.write(feature+"\tvisibility\ton\n") 
            oHandle.close()

            withClinF =1
            
    #data josn
    J={}
    J["cgDataVersion"]=1
    J["version"]= datetime.date.today().isoformat()
    J["name"]="TCGA_"+cancer+"_cohort"
    J["type"]= "clinicalMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    if withClinF:
        J[":clinicalFeature"]=  "TCGA_"+ cancer+"_cohort_clinFeat"
    oHandle = open(outfile +".json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()


    if cancer in ["LUAD","LUSC"]:
        derived_cancer="LUNG"
        doDerivedCancer(outDir, cancer, derived_cancer, flog,REALRUN)
    if cancer in ["COAD","READ"]:
        derived_cancer="COADREAD"
        doDerivedCancer(outDir, cancer, derived_cancer, flog,REALRUN)

def doDerivedCancer(outDir, cancer, derived_cancer, flog,REALRUN):
    if derived_cancer =="":
        return
    
    os.system("cp "+outDir+cancer+"/cohort "+outDir+derived_cancer+"/cohort_" +cancer)
    J={}
    J["cgDataVersion"]=1
    J["version"]= datetime.date.today().isoformat()
    J["name"]="TCGA_"+derived_cancer+"_cohort_"+cancer
    J["type"]= "clinicalMatrix" 
    J[":sampleMap"]="TCGA."+derived_cancer+".sampleMap"
    J[":clinicalFeature"]=  "TCGA_"+derived_cancer+"_cohort_"+cancer+"_clinFeat"
    oHandle = open(outDir+derived_cancer+"/cohort_" +cancer +".json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
        
    #clinFeature
    J={}
    J["name"]="TCGA_"+derived_cancer+"_cohort_"+cancer+"_clinFeat"
    J["type"]="clinicalFeature"
    oHandle = open(outDir+derived_cancer+"/cohort_" +cancer +"_clinFeature.json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
        
    oHandle = open(outDir+derived_cancer+"/cohort_" +cancer +"_clinFeature","w")
    feature="cohort"
    if TCGAUtil.featurePriority.has_key(derived_cancer):
        if TCGAUtil.featurePriority[derived_cancer].has_key(feature):
            priority= TCGAUtil.featurePriority[derived_cancer][feature]
            oHandle.write(feature+"\tpriority\t"+str(priority)+"\n")
            oHandle.write(feature+"\tvisibility\ton\n") 
    oHandle.close()
    return
