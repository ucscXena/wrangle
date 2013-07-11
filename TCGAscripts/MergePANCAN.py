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
    inFiles ={}
    outFiles={}
    for cancer in os.listdir(outDir):
        if cancer in ["LUNG","COADREAD"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename
        if not os.path.exists(cancerFile):
            continue

        #if not os.path.exists("/inside/home/jzhu/cgDataJing/scripts/TestPANCAN"+"/"+cancer):
        #    os.mkdir("/inside/home/jzhu/cgDataJing/scripts/TestPANCAN"+"/"+cancer)

        #cancerOutFile = "/inside/home/jzhu/cgDataJing/scripts/TestPANCAN"+"/"+cancer+"/"+filename+"_PANCAN"

        cancerOutFile = outDir+"/"+cancer+"/"+filename+"_PANCAN"

        inFiles[cancer]= cancerFile
        outFiles[cancer]= cancerOutFile
        
    processRNA (inFiles, filename,outFiles, REALRUN)
        
def processRNA (inFiles, filename, outFiles,REALRUN): 
    keys = inFiles.keys()

    if REALRUN:
        #header:
        for key in keys:
            fin = open(inFiles[key],'r')
            fout= open(outFiles[key],'w')
            inFiles[key] =fin
            outFiles[key]=fout
            fout.write(fin.readline())

        #data normalization per gene
        lineN=0
        while 1:
            lineN= lineN+1

            dataDic={}
            n=0
            total=0.0
            for key in keys:
                fin = inFiles[key]
                data = string.split(fin.readline(),"\t")
                for i in range(1,len(data)):
                    if data[i]=="":
                        continue
                    data[i]= float(data[i])
                    n= n+1
                    total=total + data[i]
                dataDic[key]=data

            if n ==0:
                for key in keys:
                    fout= outFiles[key]
                    fout.close()
                    outFiles[key]=fout.name
                    fin= inFiles[key]
                    fin.close()
                    inFiles[key]=fin.name
                break
        
            average = total/n

            for key in keys:
                data = dataDic[key]
                fout= outFiles[key]
                fout.write(data[0])
                for i in range(1,len(data)):
                    if data[i]=="":
                        fout.write("\t")
                        continue
                    fout.write("\t"+str(data[i]-average))
                fout.write("\n")

            print lineN
            continue
        
    for key in keys:
        cancer=key
        fin= open(inFiles[key]+".json","r")
        fout= open(outFiles[key]+".json","w")

        J= json.loads(fin.read())
        fin.close()

        J.pop("colNormalization")
        J['name']= J['name']+"_PANCAN"
        J['shortTitle']="Gene Expression (pancan normalized)"
        J["label"]= key+" Gene Expression (pancan normalized)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression (IlluminaHiSeq), pancan normalized"
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq, mean-normalized across all TCGA cohorts.<br><br>"+ \
                          " The gene expression profile was measured experimentally using the "+J['PLATFORM']+" by the "+ J['dataProducer'] +"." + \
                          " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the mean-normalized gene-level transcription estimates."
        J["wrangling_procedure"]="Level_3 Data (file names: *.rsem.genes.normalized_results) download from TCGA DCC, log2(x+1) transformed, normalized across all TCGA cancer cohorts, and deposited into UCSC cgData repository"
        J['gain'] = 0.5
        J['dataProducer']="UCSC Cancer Browser team"
        J["redistribution"]= False

        J["sample_type"]="tumor"
        J["cohort"] ="TCGA_"+cancer
        J['domain']="TCGA"
        J['owner']="TCGA"
        J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
        J["label"]= cancer +" "+J["shortTitle"]
        J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]

        fout.write(json.dumps(J,indent=-1))
        fout.close()
    return
