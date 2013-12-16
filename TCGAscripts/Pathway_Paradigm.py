import string, os, sys
import json,datetime
import inspect

LEVEL="Level_4"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmpTry/"

def Pathway_Paradigm_mRNA  (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "Pathway_Paradigm_mRNA."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm_mRNA_And_Copy_Number (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "Pathway_Paradigm_mRNA_And_Copy_Number."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)
    
def Pathway_Paradigm_RNAseq (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3]
    PATHPATTERN= "Pathway_Paradigm_RNASeq."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm_RNASeq_And_Copy_Number (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "Pathway_Paradigm_RNASeq_And_Copy_Number."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN):
    if string.find(string.upper(cancer),"PANCAN")!=-1:
        return
    
    garbage=[tmpDir]
    if os.path.exists(tmpDir ):
        os.system("rm -rf "+ tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #figure out the FH date
    if inDir[-1]!="/":
        s = inDir+"/"
    else:
        s=inDir
    FHdate = string.split(s,"/")[-2]

    #single file in dir mode, uncompress to dir
    dataDir =""
    lastDate=""
    found =0
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,"md5")==-1:
            print PATHPATTERN, file
            pass
        else:
            continue

        if not os.path.exists(inDir +file+".md5"):
            print "file has no matching .md5 throw out", file
            continue

        #file date
        lastDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)
        
        #is tar.gz?, uncompress
        if string.find(file,".tar.gz")!=-1:
            found =1
            if REALRUN:
                os.system("tar -xzf "+inDir+file +" -C "+tmpDir)
                dataDir = tmpDir +os.listdir(tmpDir)[0]+"/"
            break

    if not found:
        print "not found"
        return
    
    #make sure there is data
    if REALRUN and (dataDir =="" or (not os.path.exists(dataDir))):
    #if dataDir =="" or (REALRUN and not os.path.exists(dataDir)):
        cleanGarbage(garbage)
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    #data processing single dir mode
    foundfile=""
    pattern="inferredPathwayLevels.tab"

    if REALRUN:
        for file in os.listdir(dataDir):
            if string.find(file,pattern)!=-1:
                foundfile =file
                break

    #errorTag = cancer +"_"+__name__
    if REALRUN and foundfile =="":
        print "no file is found"
        sys.exit()

    cgFileName= PATHPATTERN[:-1]
    if REALRUN and foundfile!="":
        process(dataDir+foundfile,outDir+cancer+"/"+cgFileName,cancer)

    
    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    J={}
    #stable
    if PATHPATTERN== "Pathway_Paradigm_mRNA.":
        suffix="PDMarray"
        J["shortTitle"]= cancer +" paradigm.mRNA (Firehose)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (array mRNA)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0
        
    if PATHPATTERN== "Pathway_Paradigm_mRNA_And_Copy_Number.":
        suffix="PDMarrayCNV"
        J["shortTitle"]=cancer +" paradigm.mRNA+CNV (Firehose)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (array mRNA + CNV)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0
        
    if PATHPATTERN== "Pathway_Paradigm_RNASeq.":
        suffix="PDMRNAseq"
        J["shortTitle"]=cancer +" paradigm.RNAseq (Firehose)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (RNAseq)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0

    if PATHPATTERN=="Pathway_Paradigm_RNASeq_And_Copy_Number.":
        suffix="PDMRNAseqCNV"
        J["shortTitle"]=cancer +" paradigm.RNAseq+CNV (Firehose)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (RNAseq + CNV)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0

    J["label"] = J["shortTitle"] 
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    
    J["cgDataVersion"]=1
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= "TCGA FIREHOSE pipeline"
    J["url"]= "http://gdac.broadinstitute.org/runs/analyses__"+FHdate[0:4]+"_"+FHdate[4:6]+"_"+FHdate[6:8]+"/data/"+cancer+"/"+FHdate[0:8]+"/"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    
    #change description
    J["wrangling_procedure"]= "FIREHOSE data download from TCGA DCC, processed at UCSC into cgData repository"

    if PATHPATTERN== "Pathway_Paradigm_mRNA.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method on gene expression data from Agilent array alone."
    if PATHPATTERN== "Pathway_Paradigm_mRNA_And_Copy_Number.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method by integrating Agilent array mRNA data and copy number data."
    if PATHPATTERN=="Pathway_Paradigm_RNASeq.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method on RNAseq data alone."
    if PATHPATTERN=="Pathway_Paradigm_RNASeq_And_Copy_Number.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method by integrating RNAseq and copy number data."

    J["description"] = J["description"] +"<br><br>PARADIGM (PMID 20529912) is pathway analysis method to infer tumor sample-specific genetic activities by incorporating curated pathway interactions as well as integrating diverse types of genomic data. The pathways used in this analysis are from <a href=\"http://pid.nci.nih.gov/\" target=\"_blank\"><u>NCI pathway interaction database</u></a>."

    J["description"] = J["description"] +"<br><br>"
        
    #change cgData
    J["name"]="TCGA_"+cancer+"_"+suffix
    name = trackName_fix(J['name'])
    if name ==False:
        message = "bad object name, need fix otherwise break loader, too long "+J["name"]
        print message
        flog.write(message+"\n")
        return
    else:
        J["name"]=name
        J[":probeMap"]= "paradigmPathwayV1Probe"
        J["type"]= "genomicMatrix" 
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
                
        oHandle.write( json.dumps( J, indent=-1 ) )
        oHandle.close()

    cleanGarbage(garbage)
    return

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf "+ dir+"*")
    return

def process(infile, outfile,cancer):
    if cancer in ["LAML"]:
        cancerCode="03"
    else:
        cancerCode="01"
        
    fin=open(infile,'r')
    fout=open(outfile,'w')

    line =fin.readline()[:-1]
    ids= string.split(line,'\t')
    fout.write(ids[0])
    for id in ids[1:]:
        fout.write("\t"+id+ "-"+cancerCode)
    fout.write("\n")

    for line in fin.readlines():
        fout.write(line)
    fin.close()
    fout.close()
    
