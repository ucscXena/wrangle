import string, os, sys
import json,datetime
import csv

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from CGDataLib import *
import  TCGAUtil

def RNAseq (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return    
    print cancer, sys._getframe().f_code.co_name
    filename = "HiSeqV2"
    doAve=1
    processRNA (filename, dir,outDir, cancer,flog, REALRUN)

def Mutation (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation"
    doAve=0
    processMutation (filename, dir,outDir, cancer,flog, REALRUN)


def cohort (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "cohort"
    doAve=0
    processCohort (filename, dir,outDir, cancer,flog, REALRUN)

def processCohort (filename, dir,outDir, cancer,flog, REALRUN):
    inFiles ={}
    for cancer in os.listdir(outDir):
        if cancer in ["LUNG","COADREAD","PANCAN"]:
            continue
    
        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename
        if not os.path.exists(cancerFile):
            continue

        inFiles[cancer]= cancerFile

    cancer="PANCAN"
    outfile  = outDir+"/"+cancer+"/"+filename+"_PANCAN" 
    foutPANCAN = open(outfile,'w')
    foutPANCAN.write("sample\tcohort\n")
    
    if REALRUN:
        keys = inFiles.keys()
        for key in keys:
            fin = open(inFiles[key],'r')
            fin.readline()
            foutPANCAN.write(fin.read())
            fin.close()
        foutPANCAN.close()

    J={}
    fin= open(inFiles[keys[0]]+".json","r")
    J= json.loads(fin.read())
    fin.close()

    fout= open(outfile+".json","w")
    J["name"]= "TCGA_PANCAN_cohort"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    fout.write(json.dumps(J,indent=-1))
    fout.close()
    
def processMutation (filename, dir,outDir, cancer,flog, REALRUN):
    inFiles ={}
    outFiles={}

    for cancer in os.listdir(outDir):
        if cancer in ["LUNG","COADREAD","PANCAN"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename
        if not os.path.exists(cancerFile):
            continue

        inFiles[cancer]= cancerFile

    cancer="PANCAN"
    outFiles[cancer]= outDir+"/"+cancer+"/"+filename+"_PANCAN" 

    keys = inFiles.keys()
    
    if REALRUN:
        #header:
        foutPANCAN= open(outFiles["PANCAN"],"w")
        
        for i in range (0,len(keys)):
            fin = open(inFiles[keys[i]],'r')
            inFiles[keys[i]] =fin
            line = fin.readline()

            if i==0:
                foutPANCAN.write(string.join(string.split(line[:-1],"\t"),"\t"))
            else:
                foutPANCAN.write("\t"+string.join(string.split(line[:-1],"\t")[1:],"\t") )
        foutPANCAN.write("\n")
                
        #data normalization per gene
        lineN=0
        while 1:
            lineN= lineN+1
            end=0
            for i in range (0,len(keys)):
                fin = inFiles[keys[i]]
                line = fin.readline()
                if line =="":
                    end=1
                    inFiles[keys[i]] =fin.name
                    continue
                
                data = string.split(line[:-1],"\t")
                if i==0:
                    foutPANCAN.write(data[0])
                for i in range(1,len(data)):
                    foutPANCAN.write("\t"+data[i])

            if end:
                break

            foutPANCAN.write("\n")

            print lineN
            continue

    J={}
    fin= open(inFiles[keys[0]]+".json","r")
    J= json.loads(fin.read())
    fin.close()

    cancer="PANCAN"
    fout= open(outFiles[cancer]+".json","w")
    
    if filename=="mutation":
        mutationJSON(J,cancer)

    J['dataProducer']="UCSC Cancer Browser team"
    J["redistribution"]= False
    J["sample_type"]="tumor"
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['owner']="TCGA"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    
    J["primary_disease"]="cancer"
    J["anatomical_origin"]= ""
        
    fout.write(json.dumps(J,indent=-1))
    fout.close()
    return

def mutationJSON(J, cancer):
    J['name']= "TCGA_PANCAN_mutation"
    J["label"]= cancer +" mutation"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" somatic mutation"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation. Red color (=1) represents non-silent somatic mutations (nonsense, missense, frame-shift indels, splice site mutations, stop codon readthroughs) are identified in the protein coding region of a gene, and white color (=0) means that none of the above mutation calls are made in this gene for the specific sample."
    J["wrangling_procedure"]="Data is combined from all TCGA cohorts and deposited into UCSC cgData repository"
    return
    
def processRNA (filename, dir,outDir, cancer,flog, REALRUN):
    inFiles ={}
    outFiles={}
    for cancer in os.listdir(outDir):
        if cancer in ["LUNG","COADREAD","PANCAN"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename
        if not os.path.exists(cancerFile):
            continue

        cancerOutFile = outDir+"/"+cancer+"/"+filename+"_PANCAN"

        inFiles[cancer]= cancerFile
        outFiles[cancer]= cancerOutFile

    cancer="PANCAN"
    outFiles[cancer]= outDir+"/"+cancer+"/"+filename+"_PANCAN" 
        
    keys = inFiles.keys()

    if REALRUN:
        #header:
        foutPANCAN= open(outFiles["PANCAN"],"w")

        for i in range (0,len(keys)):
            fin = open(inFiles[keys[i]],'r')
            fout= open(outFiles[keys[i]],'w')
            inFiles[keys[i]] =fin
            outFiles[keys[i]]=fout
            line = fin.readline()
            fout.write(line)

            if i==0:
                foutPANCAN.write(string.join(string.split(line[:-1],"\t"),"\t"))
            else:
                foutPANCAN.write("\t"+string.join(string.split(line[:-1],"\t")[1:],"\t") )
        foutPANCAN.write("\n")
                
        #data normalization per gene
        lineN=0
        while 1:
            lineN= lineN+1

            dataDic={}

            n=0
            total=0.0
            end=0
            for key in keys:
                fin = inFiles[key]
                line = fin.readline()
                if line =="":
                    end =1
                    inFiles[key] = fin.name
                    outFiles[key] = fout.name
                    continue
                
                data = string.split(line[:-1],"\t")
                
                for i in range(1,len(data)):
                    if data[i]=="":
                        continue
                    data[i]= float(data[i])
                    n= n+1
                    total=total + data[i]
                dataDic[key]=data
            if end:
                break
            if n==0:
                average=0
            else:
                average = total/n

            for i in range (0,len(keys)):
                data = dataDic[keys[i]]
                fout= outFiles[keys[i]]
                fout.write(data[0])
                if i==0:
                    foutPANCAN.write(data[0])
                for i in range(1,len(data)):
                    if data[i]=="":
                        fout.write("\t")
                        foutPANCAN.write("\t")
                        continue
                    fout.write("\t"+str(data[i]-average))
                    foutPANCAN.write("\t"+str(data[i]-average))
                fout.write("\n")

            foutPANCAN.write("\n")

            print lineN
            continue

    keys.append("PANCAN")
    J={}
    for key in keys:
        cancer=key

        if inFiles.has_key(key):
            fin= open(inFiles[key]+".json","r")
            J= json.loads(fin.read())
            fin.close()

        fout= open(outFiles[key]+".json","w")

        if J.has_key("colNormalization"):
            J.pop("colNormalization")

        J['name']= J['name']+"_PANCAN"
        J['dataProducer']="UCSC Cancer Browser team"
        J["redistribution"]= False
        J["sample_type"]="tumor"
        J["cohort"] ="TCGA_"+cancer
        J['domain']="TCGA"
        J['owner']="TCGA"
        J["label"]= cancer +" "+J["shortTitle"]
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
        J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]
        
        J['shortTitle']="Gene Expression (pancan normalized)"
        J["label"]= key+" gene expression (pancan normalized)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq (IlluminaHiSeq), pancan normalized"
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq, mean-normalized across all TCGA cohorts.<br><br>"+ \
                          " The gene expression profile was measured experimentally using the "+J['PLATFORM']+" by the "+ J['dataProducer'] +"." + \
                          " This dataset shows the mean-normalized gene-level transcription estimates."
        J["wrangling_procedure"]="Level_3 Data (file names: *.rsem.genes.normalized_results) download from TCGA DCC, log2(x+1) transformed, normalized across all TCGA cancer cohorts, and deposited into UCSC cgData repository"
        J['gain'] = 0.5
        
        if cancer!="PANCAN":
            J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
            J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
        else:
            J["primary_disease"]="cancer"
            J["anatomical_origin"]= ""
            J["url"]="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
            J["name"]= "TCGA_PANCAN_exp_HiSeqV2_PANCAN"
        fout.write(json.dumps(J,indent=-1))
        fout.close()
    return
