import sys,string,os
import json,datetime
import math
import inspect
import copy

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmptmp/"

#/inside/home/cline/projects/PanCancer/mutationMatrices/*_cleaned_filtered.txt

def mutationMatrix (inDir, outDir, cancer,flog,REALRUN):
    if string.find(cancer,"PANCAN") !=-1:
        return
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "_cleaned_filtered.txt"
    namesuffix = "mutation"
    dataProducer = "TCGA PANCAN"

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
    J[":dataSubType"]="mutation"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]="https://www.synapse.org/#!Synapse:syn1729383"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J[":probeMap"]= "hugo"
    J["shortTitle"]="Mutation"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation"
    J["gain"]=10
    J["description"]= "The dataset shows TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data. Red color (=1) represents non-silent somatic mutations (nonsense, missense, frame-shift indels, splice site mutations, stop codon readthroughs) are identified in the protein coding region of a gene, and white color (=0) means that none of the above mutation calls are made in this gene for the specific sample. Somatic mutations calls (even on the same tumor DNA extract) are effected by many factors including library prep, sequencing process, reads mapping method, reference genome used, calling algorithms, and ad-hoc pre/postprocessing such as black list genes, target selection regions, and black list samples.  This dataset is the best effort made by the TCGA PANCANER analysis working group.  PANCAN mutation data can be downloade at the url site shown below."
    J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
    J["wrangling_procedure"]= "TCGA PANCAN strictly filtered maf files (file names: *_cleaned_filtered.maf) download from Synapse, processed into gene by sample matrix at UCSC into cgData repository"

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

