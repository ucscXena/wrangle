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


def clin (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return

    print cancer, sys._getframe().f_code.co_name
    filename = "_clinicalMatrix"
    doAve=0
    dir = "/inside/home/jzhu/cgDataJing/scripts/data_flatten/public/TCGA/"
    processClin (filename, dir,outDir, cancer,flog, REALRUN)

def processClin (filename, dir,outDir, CANCER,flog, REALRUN):
    inFiles ={}
    print os.path.dirname(dir)
    for cancer in os.listdir(os.path.dirname(dir)):
        if cancer in ["LUNG","COADREAD","PANCAN"]:
            continue

        cancerDir= os.path.dirname(dir)  + "/"+cancer
        cancerFile = cancerDir+"/"+cancer+filename
        if not os.path.exists(cancerFile):
            continue

        inFiles[cancer]= cancerFile

    features=["sample_type","sample_type_id","gender","cohort","age_at_initial_pathologic_diagnosis",
              "_OS","_OS_IND","_RFS","_RFS_IND","_EVENT","_TIME_TO_EVENT"]
        
    for feature in features:
        print feature
        if REALRUN:
            outfile  = outDir+"/"+CANCER+"/"+feature+"_PANCAN"
            foutPANCAN = open(outfile,'w')
            foutPANCAN.write("sample\t"+feature+"\n")
            
            keys = inFiles.keys()
            os.system("rm -f tmp")
            samples=[]

            for key in keys:
                file = inFiles[key]
                fin = open(file,'r')
                POS = 0
                data = string.split(fin.readline()[:-1],"\t")
                for i in range (0,len(data)):
                    if data[i] == feature:
                       POS=i
                       break
                if POS==0:
                    continue
                #print file, POS
                for line in fin.readlines():
                    data = string.split(line[:-1],"\t")
                    sample=data[-1]
                    if sample not in samples and sample!="":
                        samples.append(sample)
                        foutPANCAN.write(data[-1]+"\t"+data[POS]+"\n")
            foutPANCAN.close()

        outfile  = outDir+"/"+CANCER+"/"+feature+"_PANCAN.json"
        fout =open(outfile,'w')

        J={}
        J['name']=feature+"_PANCAN"
        J['type']="clinicalMatrix"
        J[":sampleMap"]="TCGA."+CANCER+".sampleMap"

        if TCGAUtil.featurePriority.has_key(CANCER) or TCGAUtil.valueType.has_key(feature):
            featureConfig=0
            cfout =open(outDir+"/"+CANCER+"/"+feature+"_PANCAN_clinFeature","w")
            
            if TCGAUtil.featurePriority[CANCER].has_key(feature):
                featureConfig=1
                priority= TCGAUtil.featurePriority[CANCER][feature]
                cfout.write(feature+"\tpriority\t"+str(priority)+"\n")
                cfout.write(feature+"\tvisibility\ton\n")

            stateOrder=None
            if TCGAUtil.featureStateOrder.has_key(feature):
                if TCGAUtil.featureStateOrder[feature].has_key(CANCER):
                    featureConfig=1
                    stateOrder = TCGAUtil.featureStateOrder[feature][CANCER]
                if TCGAUtil.featureStateOrder[feature].has_key("ALL"):
                    featureConfig=1
                    stateOrder = TCGAUtil.featureStateOrder[feature]["ALL"]
            if stateOrder:
                cfout.write(feature+"\tvalueType\tcategory\n")
                for state in stateOrder:
                    cfout.write(feature+"\tstate\t"+state+"\n")

                cfout.write(feature+"\tstateOrder\t\""+string.join(stateOrder,"\",\"")+"\"\n")
                cfout.write(feature+"\tstateOrderRelax\ttrue\n")

            if TCGAUtil.valueType.has_key(feature):
                featureConfig=1
                cfout.write(feature+"\tvalueType\tcategory\n")
                
            cfout.close()
            
            if  featureConfig:
                cfJ ={}
                cfJ["name"]= J['name']+"_clinFeature"
                cfJ["type"]="clinicalFeature"
                cfout=open(outDir+"/"+CANCER+"/"+feature+"_PANCAN_clinFeature.json","w")
                cfout.write(json.dumps(cfJ,indent=-1))
                cfout.close()
                
                J[":clinicalFeature"]= cfJ["name"]
                
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
    J["sample_type"]=["tumor"]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"] 
    J['owner']="TCGA"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]

            
    if cancer!="PANCAN":
        J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
        J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    else:
        J["primary_disease"]="cancer"
        J["anatomical_origin"]=""
        origin =[]
        for value in TCGAUtil.anatomical_origin.values():
            if value =="":
                continue
            if value not in origin:
                origin.append(value)
        J["tags"]= J["tags"]+origin
        
    fout.write(json.dumps(J,indent=-1))
    fout.close()
    return

def mutationJSON(J, cancer):
    J['name']= "TCGA_PANCAN_mutation"
    J["shortTitle"]= cancer +" mutation"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" somatic mutation"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation. Red color (=1) represents non-silent somatic mutations (nonsense, missense, frame-shift indels, splice site mutations, stop codon readthroughs) are identified in the protein coding region of a gene, and white color (=0) means that none of the above mutation calls are made in this gene for the specific sample."
    J["description"] = J["description"] +"<br><br>"
    
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
                    fout =outFiles[key]
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
        J["sample_type"]="tumor"
        J["cohort"] ="TCGA_"+cancer
        J['domain']="TCGA"
        J['owner']="TCGA"
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
        J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]
        
        J['shortTitle']=key+" gene expression (pancan normalized)"
        J["label"] = J["shortTitle"] 
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
            J["anatomical_origin"]=""
            origin =[]
            for value in TCGAUtil.anatomical_origin.values():
                if value =="":
                    continue
                if value not in origin:
                    origin.append(value)

            J["tags"]= J["tags"]+origin
            J["url"]="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
            J["name"]= "TCGA_PANCAN_exp_HiSeqV2_PANCAN"
        fout.write(json.dumps(J,indent=-1))
        fout.close()
    return
