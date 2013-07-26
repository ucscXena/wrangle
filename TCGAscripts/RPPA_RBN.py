import sys,string,os
import json,datetime
import math
import inspect
import copy

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmptmp/"

#/inside/home/cline/projects/PanCancer/mutationMatrices/*_RBN

def RPPA_RBN (inDir, outDir, cancer,flog,REALRUN):
    if string.find(cancer,"PANCAN") !=-1:
        return
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "_RBN"
    namesuffix = "RPPA_RBN"
    dataProducer = "MDACC"

    garbage=[tmpDir]

    if os.path.exists( tmpDir ):
        os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #inDir is a single file
    rootDir =""
    remoteDataDirExample =""
    #find the file
    file = inDir

    if string.find(file,PATHPATTERN)!=-1 and string.find(string.upper(file),cancer)!=-1: 
        pass
    else:
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= namesuffix 

    if REALRUN:
        os.system("cp "+file +" "+outDir+cancer+"/"+cgFileName)

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")

    J={}
    #stable
    J["cgDataVersion"]=1
    J[":dataSubType"]="RPPA"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]="https://www.synapse.org/#!Synapse:syn1759392"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J[":probeMap"]= "md_anderson_antibodies"
    J["shortTitle"]="Protein (RBN)"
    J["label"]= cancer +" protein (RBN)"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") reverse phase protein array (replicate-base normalization)"

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]="tumor"
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['owner']="TCGA"
    
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") protein expression data for 131 proteins, measured by RPPA (reverse phase protein array) technology. These data have been normalized by RBN (replicate-base normalization) method developed by MDACC. Details: https://www.synapse.org/#!Synapse:syn1750330.<br><br>"

    J["wrangling_procedure"]= "Data download from https://www.synapse.org/#!Synapse:syn1759392, processed into antibody by sample matrix at UCSC into cgData repository"

    #change cgData
    J["name"]="TCGA_"+cancer+"_"+namesuffix
    name = trackName_fix(J['name'])
    if name ==False:
        message = "bad object name, need fix otherwise break loader, too long "+J["name"]
        print message
        flog.write(message+"\n")
        return
    else:
        J["name"]=name        
        
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
    
    return

