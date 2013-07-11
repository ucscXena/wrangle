import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def cohort (inDir, outDir, cancer, flog,REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name
    if cancer in ["COADREAD","LUNG"]:
        return
        
    file =inDir+"/"+cancer+"_clinicalMatrix"
    outfile = outDir+cancer+"/cohort"
    os.system("echo \"sample\tcohort\" > "+outfile)
    os.system("c=$(head -n 1 "+file +"| tr \"\t\" \"\n\" | wc -l); cut -f $c "+file +" | sort |uniq | grep -v _INTEGRATION | grep -v ^$|awk 'BEGIN{FS=OFS=\"\t\"; cohort=\""+cancer+"\"} {print $1,cohort}' >> "+outfile )

    J={}
    J["cgDataVersion"]=1
    J["version"]= datetime.date.today().isoformat()
    J["name"]="TCGA_"+cancer+"_cohort"
    J["type"]= "clinicalMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    oHandle = open(outfile +".json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
    
    derived_cancer =""
    if cancer in ["LUAD","LUSC"]:
        derived_cancer="LUNG"
    if cancer in ["COAD","READ"]:
        derived_cancer="COADREAD"

    if derived_cancer !="":
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
