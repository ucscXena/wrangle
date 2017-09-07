import sys,string,os
import json,datetime
import math
import inspect
import copy

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmpTry/"

#/inside/home/cline/projects/PanCancer/mutationMatrices/*_RBN

def RPPA_RBN (inDir, outDir, cancer,flog,REALRUN):
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
    J["dataSubType"]="protein expression RPPA"
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]="https://www.synapse.org/#!Synapse:syn1759392"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J[":probeMap"]= "md_anderson_antibodies"
    J["label"]= "RPPA (replicate-base normalization)"
    J["sample_type"]=["tumor"]
    J["redistribution"]= True
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    if cancer !="PANCAN12":
        J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
        J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") phospho- or total protein expression by reverse phase protein array (replicate-base normalization)"
    else:
        J["primary_disease"] = "cancer"
        J["anatomical_origin"]= ""
        J["longTitle"]="TCGA PANCAN AWG phospho- or total- protein expression by reverse phase protein array (replicate-base normalization)"
        J["citation"] = "Cell 2014 http://dx.doi.org/10.1016/j.cell.2014.06.049"
        J["articleTitle"]= "Multi-platform analysis of 12 cancer types reveals molecular classification within and across tissues-of-origin"
        J["dataProducer"]= "Hoadley, K.M., Yau, C., Wolf, D.M., et al. and TCGA PANCAN AWG"

        J["tags"]=["PANCAN12","PANCAN11"]
        origin =[]
        for values in [ TCGAUtil.anatomical_origin["BLCA"],["BLCA"],\
                            TCGAUtil.anatomical_origin["BRCA"],["BRCA"],\
                            TCGAUtil.anatomical_origin["COAD"],["COAD"],\
                            TCGAUtil.anatomical_origin["UCEC"],["UCEC"],\
                            TCGAUtil.anatomical_origin["GBM"],["GBM"],\
                            TCGAUtil.anatomical_origin["HNSC"],["HNSC"],\
                            TCGAUtil.anatomical_origin["KIRC"],["KIRC"],\
                            TCGAUtil.anatomical_origin["LUAD"],["LUAD"],\
                            TCGAUtil.anatomical_origin["LUSC"],["LUSC"],\
                            TCGAUtil.anatomical_origin["OV"],["OV"],\
                            TCGAUtil.anatomical_origin["READ"],["READ"]
                        ]:
            for value in values:
                if value not in origin:
                    origin.append(value)

        J["tags"]= J["tags"]+origin
      
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['owner']="TCGA"
    
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") protein expression data for 131 proteins, measured by RPPA (reverse phase protein array) technology. These data have been normalized by RBN (replicate-base normalization) method developed by MDACC. Details: https://www.synapse.org/#!Synapse:syn1750330 and http://bioinformatics.mdanderson.org/main/TCPA:Overview."
    J["description"] = J["description"] +"<br><br>"
    J["unit"]="RBN normalized RPPA value"
    J["wrangling_procedure"]= "Data download from https://www.synapse.org/#!Synapse:syn1759392, processed into antibody by sample matrix at UCSC into Xena repository"

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

