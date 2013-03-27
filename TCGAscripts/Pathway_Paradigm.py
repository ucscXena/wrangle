import string, os, sys
import json,datetime
import inspect

LEVEL="Level_4"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

def Pathway_Paradigm_agilent (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "ParadigmReport."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm_copynumber (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "ParadigmReportWithCopyNumber."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)
    
def Pathway_Paradigm_RNAseq (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3]
    PATHPATTERN= "ParadigmReportWithRNASeq."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm_RNAseq_copynumber (inDir, outDir, cancer,flog, REALRUN):
    #print status
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "ParadigmReportWithRNASeqAndCopyNumber."
    Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN)

def Pathway_Paradigm (inDir, outDir, cancer,flog, PATHPATTERN, REALRUN):
    garbage=["tmptmp/"]
    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")

    #figure out the FH date
    if inDir[-1]!="/":
        s = inDir+"/"
    else:
        s=inDir
    FHdate = string.split(s,"/")[-2]

    #single file in dir mode, uncompress to dir
    dataDir =""
    lastDate=""
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
            if REALRUN:
                os.system("tar -xzf "+inDir+file +" -C tmptmp/")
                dataDir = "tmptmp/"+os.listdir("tmptmp/")[0]+"/"
            break
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
    if PATHPATTERN== "ParadigmReport.":
        suffix="PDMarray"
        J["shortTitle"]="Paradigm.array_expression"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (agilent expression)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0
        
    if PATHPATTERN== "ParadigmReportWithCopyNumber.":
        suffix="PDMCNV"
        J["shortTitle"]="Paradigm.CNV"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (CNV)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0
        
    if PATHPATTERN== "ParadigmReportWithRNASeq.":
        suffix="PDMExp"
        J["shortTitle"]="Paradigm.RNAseq"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (RNAseq)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0

    if PATHPATTERN=="ParadigmReportWithRNASeqAndCopyNumber.":
        suffix="PDMExpCNV"
        J["shortTitle"]="Paradigm.RNAseq+CNV"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") PARADIGM inference activity (RNAseq,CNV)"
        J[":dataSubType"]="PARADIGM"
        J["gain"]=1.0

                        
    J["cgDataVersion"]=1
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= "TCGA FIREHOSE pipeline"
    J["url"]= "http://gdac.broadinstitute.org/runs/analyses__"+FHdate[0:4]+"_"+FHdate[4:6]+"_"+FHdate[6:8]+"/data/"+cancer+"/"+FHdate[0:8]+"/gdac.broadinstitute.org_"+cancer+"."+PATHPATTERN+"Level_4."+FHdate+".0.0.tar.gz"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    
    #change description
    J["wrangling_procedure"]= "FIREHOSE data download from TCGA DCC, processed at UCSC into cgData repository"

    if PATHPATTERN== "ParadigmReport.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method on gene expression data from Agilent array alone."
    if PATHPATTERN== "ParadigmReportWithCopyNumber.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method on copy number data alone."
    if PATHPATTERN=="ParadigmReportWithRNASeq.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method on RNAseq data alone."
    if PATHPATTERN=="ParadigmReportWithRNASeqAndCopyNumber.":
        J["description"]= "Broad FireHose automated run results of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"+\
                          " gene activity level inferred using the PARADIGM method by integrating RNAseq and copy number data."

    J["description"] = J["description"] +" PARADIGM is pathway analysis method to infer patient or sample-specific genetic activities by incorporating curated pathway interactions as well as integrating diverse types of genomic data, such as gene expression and copy number data. The pathways used in this analysis are from <a href=\"http://pid.nci.nih.gov/\" target=\"_blank\"><u>NCI pathway interaction database</u></a>."+\
                       " Genes are mapped onto the human genome coordinates using cgData HUGO probeMap."+\
                       " Reference to PARADIGM pathway analysis method: PMID 20529912."
    J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
        
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
        os.system("rm -rf dir")
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
    
