import string, os, sys
import json,datetime
import csv

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *
from CGDataUtil import *
from CGDataLib import *
import  TCGAUtil
import mergeGenomicMatrixFiles
import xenaToMatrix

def RNAseq (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return    
    print cancer, sys._getframe().f_code.co_name
    filename = "HiSeqV2"
    processRNA(filename, dir,outDir, cancer,flog, REALRUN)

def HiSeqV2  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "HiSeqV2"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

def HiSeqV2_exon  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "HiSeqV2_exon"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

def miRNA  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "miRNA"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

def mutation_broad  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_broad"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def mutation_wustl  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_wustl"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def mutation_bcgsc  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_bcgsc"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def mutation_bcm  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_bcm"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def mutation_ucsc_vcf  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_ucsc_vcf"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def Gistic2 (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "Gistic2_CopyNumber_Gistic2_all_data_by_genes"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

def clin (dir,outDir, cancer,flog,REALRUN):
    if cancer not in ["PANCAN12","PANCAN"]:
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "_clinical"
    doAve=0
    dir = "/inside/home/jzhu/cgDataJing/scripts/data_flatten/public/TCGA/"
    processClin (filename, dir,outDir, cancer,flog, REALRUN)
    
def processClin (filename, dir,outDir, CANCER,flog, REALRUN):
    inFiles ={}
    for cancer in os.listdir(os.path.dirname(dir)):
        if cancer in ["LUNG","COADREAD","GBMLGG", "PANCAN","PANCAN12"]:
            continue

        cancerDir= os.path.dirname(dir)  + "/"+cancer

        cancerFile = cancerDir+"/"+cancer+filename+"Matrix"

        if not os.path.exists(cancerFile):
            continue

        inFiles[cancer]= cancerFile

    features=["sample_type","sample_type_id","gender","_cohort","age_at_initial_pathologic_diagnosis",
              "_OS","_OS_IND","_OS_UNIT",
              "_RFS","_RFS_IND","_RFS_UNIT",
              "_EVENT","_TIME_TO_EVENT","_TIME_TO_EVENT_UNIT",
              "_anatomical_origin","_primary_disease"]   


    for feature in features:
        print feature
        if REALRUN:
            outfile  = outDir+"/"+CANCER+"/"+feature+"_"+CANCER
            foutPANCAN = open(outfile,'w')
            foutPANCAN.write("sample\t"+feature+"\n")
            
            keys = inFiles.keys()
            os.system("rm -f tmp")
            samples=[]

            for key in keys:
                file = inFiles[key]

                # check if it is deprecated
                featureFile = os.path.dirname(dir)  + "/"+key+"/"+key+filename+"Feature"
                if os.path.exists(featureFile):
                    clinFeature = ClinicalFeatureNew.ClinicalFeatureNew(featureFile,"tmpName")
                    longTitle= clinFeature.getLongTitle(feature)
                    if longTitle and string.find(longTitle,"_DEPRECATED_")!=-1:
                        continue
                
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

        outfile  = outDir+"/"+CANCER+"/"+feature+"_"+CANCER+".json"
        fout =open(outfile,'w')

        J={}
        J['name']=feature+"_"+CANCER
        J['type']="clinicalMatrix"
        J[":sampleMap"]="TCGA."+CANCER+".sampleMap"
        J["dataSubType"]="phenotype"

        if TCGAUtil.featurePriority.has_key(CANCER) or TCGAUtil.valueType.has_key(feature):
            featureConfig=0
            
            if TCGAUtil.featurePriority[CANCER].has_key(feature):
                cfout =open(outDir+"/"+CANCER+"/"+feature+"_"+CANCER+"_clinFeature","w")
                cfout.write("#feature\tattribute\tvalue\n")
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
                if featureConfig==0:
                    cfout =open(outDir+"/"+CANCER+"/"+feature+"_"+CANCER+"_clinFeature","w")
                    cfout.write("#feature\tattribute\tvalue\n")
                featureConfig=1
                cfout.write(feature+"\tvalueType\tcategory\n")
                for state in stateOrder:
                    cfout.write(feature+"\tstate\t"+state+"\n")

                cfout.write(feature+"\tstateOrder\t\""+string.join(stateOrder,"\",\"")+"\"\n")
                cfout.write(feature+"\tstateOrderRelax\ttrue\n")

            if TCGAUtil.valueType.has_key(feature):
                if featureConfig==0:
                    cfout =open(outDir+"/"+CANCER+"/"+feature+"_"+CANCER+"_clinFeature","w")
                    cfout.write("#feature\tattribute\tvalue\n")
                featureConfig=1
                cfout.write(feature+"\tvalueType\tcategory\n")
                
            cfout.close()
            
            if  featureConfig:
                cfJ ={}
                cfJ["name"]= J['name']+"_clinFeature"
                cfJ["type"]="clinicalFeature"
                cfout=open(outDir+"/"+CANCER+"/"+feature+"_"+CANCER+"_clinFeature.json","w")
                cfout.write(json.dumps(cfJ,indent=-1))
                cfout.close()
                
                J[":clinicalFeature"]= cfJ["name"]
                
        fout.write(json.dumps(J,indent=-1))
        fout.close()
        
def processFiles (filename, dir,outDir, cancer , MUTATION= False):
    inFiles =[]

    for cancer in os.listdir(outDir):
        if cancer in ["LUNG","COADREAD","PANCAN","GBMLGG","PANCAN12"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename

        if not os.path.exists(cancerFile):
            if filename =="miRNA":
                cancerFile = cancerDir+"/"+filename+"_HiSeq"
                if not os.path.exists(cancerFile):
                    cancerFile = cancerDir+"/"+filename+"_GA"
                    if not os.path.exists(cancerFile):
                        continue
                    else:
                        inFiles.append(cancerFile)
                else:
                    inFiles.append(cancerFile)
            else:
                continue

        if MUTATION:
            fin = open(cancerFile+".json",'r')
            J= json.loads(fin.read())
            if (not J.has_key("assembly")) or J["assembly"]!="hg19":
                continue

        inFiles.append(cancerFile)
    return inFiles

def convertToMatrix(infile,cancer):
    samples=[]
    genes=[]
    dic={}
    outfile = infile+"_gene"
    allGeneList = xenaToMatrix.allGeneList(xenaToMatrix.ALLGENE_FILE)
    xenaToMatrix.process (infile, samples, allGeneList, dic) 
    xenaToMatrix.output(outfile, samples, genes, dic)

def processMutationData (filename, dir, outDir, cancer,flog, REALRUN):
    inFiles =processFiles (filename, dir, outDir, cancer , True)

    cancer="PANCAN"
    outfile = outDir +cancer + "/"+filename

    if REALRUN:
        fout= open(outfile,"w")
        header =""
        for infile in inFiles:
            fin = open(infile,'r')
            if header =="":
                header = fin.readline()
                fout.write(header)
            else:
                if header != fin.readline():
                    print "check header"
            fout.write(fin.read())
                    
            fin.close()
        fout.close()

        convertToMatrix(outfile,cancer)

    J={}
    fout = open(outfile+".json","w")

    cancer="PANCAN"    
    J["type"]= "mutationVector"
    genelevel = ""
    if filename=="mutation_broad":
        mutation_broadJSON(J,cancer,genelevel)
    elif filename=="mutation_bcm":
        mutation_bcmJSON(J,cancer,genelevel)
    elif filename=="mutation_ucsc_vcf":
        mutation_ucsc_vcfJSON(J,cancer,genelevel)
    elif filename=="mutation_bcgsc":
        mutation_bcgscJSON(J,cancer,genelevel)
    elif filename=="mutation_wustl":
        mutation_wustlSON(J,cancer,genelevel)
    commonJSON(J, cancer)
    J["dataSubType"]="somatic mutation (SNPs and small INDELs)"
    fout.write(json.dumps(J,indent=-1))
    fout.close()


    J={}
    fout = open(outfile+"_gene.json","w")

    cancer="PANCAN"    
    J["type"]= "genomicMatrix"
    genelevel = " gene-level"
    if filename=="mutation_broad":
        mutation_broadJSON(J,cancer,genelevel)
    elif filename=="mutation_bcm":
        mutation_bcmJSON(J,cancer,genelevel)
    elif filename=="mutation_ucsc_vcf":
        mutation_ucsc_vcfJSON(J,cancer,genelevel)
    elif filename=="mutation_bcgsc":
        mutation_bcgscJSON(J,cancer,genelevel)
    elif filename=="mutation_wustl":
        mutation_wustlSON(J,cancer,genelevel)
    J[":probeMap"]="hugo"
    J["dataSubType"]="somatic non-silent mutation (gene-level)"
    commonJSON(J, cancer)

    fout.write(json.dumps(J,indent=-1))
    fout.close()

    return

def processMatrix (filename, dir,outDir, cancer,flog, REALRUN):
    inFiles =processFiles (filename, dir,outDir, cancer )

    cancer="PANCAN"
    outfile = outDir +cancer + "/"+filename

    if REALRUN:
        #header:
        foutPANCAN= open(outfile,"w")
        
        genes={}
        samples={}
        dataMatrix=[]

        for infile in inFiles:
            mergeGenomicMatrixFiles.header (samples, infile)
        for infile in inFiles:
            mergeGenomicMatrixFiles.process(genes, samples, dataMatrix, infile)
        mergeGenomicMatrixFiles.outputMatrix(dataMatrix, samples, genes, outfile)                

    J={}
    fout = open(outfile+".json","w")

    cancer="PANCAN"    
    J["type"]= "genomicMatrix"
    if filename=="Gistic2_CopyNumber_Gistic2_all_data_by_genes":
        gisticJSON(J,cancer)
    elif filename=="HiSeqV2":
        HiSeqV2JSON(J,cancer)
    elif filename=="HiSeqV2_exon":
        HiSeqV2_exonJSON(J,cancer)
    elif filename=="miRNA":
        miRNAJSON(J,cancer)
 
    commonJSON(J, cancer)

    fout.write(json.dumps(J,indent=-1))
    fout.close()
    return

def commonJSON(J, cancer):
    J['dataProducer']="UCSC Xena team"
    J["sample_type"]=["tumor"]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["wrangling_procedure"]="Data is combined from all TCGA cohorts and deposited into UCSC Xena repository"
    J["version"]= datetime.date.today().isoformat()
    J["primary_disease"]="cancer"
    J["anatomical_origin"]=""
    origin =[]
    for value in TCGAUtil.anatomical_origin.values():
        for data in value:
            if data not in origin:
                origin.append(data)
    J["tags"]= J["tags"]+origin
    return


def mutation_bcgscJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_bcgsc"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (bcgsc)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" somatic mutation (bcgsc)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (bcgsc)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" nonsilent somatic mutation (bcgsc)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data. The calls are generated at Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC) using the BCGSC pipeline method. BCGSC's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    

def mutation_wustlSON (J, cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_wustl"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (wustl)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" somatic mutation (wustl)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (wustl)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" nonsilent somatic mutation (wustl)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data. The calls are generated at Washington University Genome Center using the WashU pipeline method. BCGSC's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"

def mutation_ucsc_vcfJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_ucsc_vcf"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (ucsc automated vcf)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" somatic mutation (ucsc)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (ucsc automated vcf)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" nonsilent somatic mutation (ucsc)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data.  The calls are generated at University of Californis Santa Cruz GDAC using the RADIA method. RADIA's calls (vcfs) from various TCGA cohorts are combined to produce this dataset. Reference to RADIA: PMID: 25405470."
    J["description"] = J["description"] +"<br><br>"    

def mutation_bcmJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_bcm"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (bcm)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+ " somatic mutation (bcm)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (bcm)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+ " nonsilent somatic mutation (bcm)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data.  The calls are generated at Baylor College of Medicine Human Genome Sequencing Center using the Baylor pipeline method. Baylor's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    

def mutation_broadJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_broad"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (broad)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" somatic mutation (broad)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (broad)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" nonsilent somatic mutation (broad)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data.  The calls are generated at Broad Institute Genome Sequencing Center using the MuTect method. MuTect calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    
    
def miRNAJSON(J,cancer):
    J['name']= "TCGA_PANCAN_miRNA"
    J["shortTitle"]= "miRNA expression"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" miRNA expression (RNAseq)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" miRNA expression by RNAseq.<br>"
    J["description"] = J["description"] +"miRNA expression measured using the IlluminaHiSeq and IllunimaGA technology. Data from all TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]= "miRBase_primary_transcript_hg18"
    J["dataSubType"]="miRNA expression RNAseq"
    J["colNormalization"]=True
    return

def HiSeqV2JSON (J, cancer):
    J['name']= "HiSeqV2_PANCAN"
    J["shortTitle"]= "gene expression"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" gene expression (IlluminaHiSeq)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" gene expression by RNAseq.<br>"
    J["description"] = J["description"] +"Gene expression measured using the IlluminaHiSeq technology. Data from all TCGA cohorts are combined to produce this dataset. Values are log2(x+1) transformed RSEM gene-level expression estimations."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="hugo"
    J["dataSubType"]="gene expression RNAseq"
    J["colNormalization"]=True
    J["wrangling_procedure"] = "Level_3 data (file names: *.rsem.genes.normalized_results) are downloaded from each cancer project at TCGA DCC, log2(x+1) transformed, and then combined at UCSC into Xena repository."
    return

def HiSeqV2_exonJSON (J, cancer):
    J['name']= "HiSeqV2_exon_PANCAN"
    J["shortTitle"]= "exon expression"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" exon expression (IlluminaHiSeq)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" exon expression by RNAseq.<br>"
    J["description"] = J["description"] +"Exon expression measured using the IlluminaHiSeq technology. Data from all TCGA cohorts are combined to produce this dataset. Values are log2(x+1) transformed exon-level transcription estimates in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="unc_RNAseq_exon"
    J["dataSubType"]="exon expression RNAseq"
    J["colNormalization"]=True
    J["wrangling_procedure"] = "Level_3 data (file names: *.exon_quantification.txt) are downloaded from each cancer project at TCGA DCC, log2(x+1) transformed, and then combined at UCSC into Xena repository."
    return


def gisticJSON(J,cancer):
    J['name']= "TCGA_PANCAN_gistic2"
    J["shortTitle"]= "gene-level copy number (gistic2)"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" gene-level copy number (gistic2)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" gene-level copy number variation (CNV) estimated using the GISTIC2 method.<br>"
    J["description"] = J["description"] +"Copy number profile was measured experimentally using whole genome microarray at Broad TCGA genome characterization center. Subsequently, TCGA FIREHOSE pipeline applied GISTIC2 method to produce segmented CNV data, which was then mapped to genes to produce gene-level estimates. Gistic2 data from all TCGA cohorts are combined to produce this dataset. Reference to GISTIC2 method PMID:21527027."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="hugo"
    J["dataSubType"]="copy number"
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

    #cancer="PANCAN"
    #outFiles[cancer]= outDir+"/"+cancer+"/"+filename+"_PANCAN" 
        
    keys = inFiles.keys()

    if REALRUN:
        #header:
        #foutPANCAN= open(outFiles["PANCAN"],"w")

        for i in range (0,len(keys)):
            fin = open(inFiles[keys[i]],'r')
            fout= open(outFiles[keys[i]],'w')
            inFiles[keys[i]] =fin
            outFiles[keys[i]]=fout
            line = fin.readline()
            fout.write(line)

            #if i==0:
            #    foutPANCAN.write(string.join(string.split(line[:-1],"\t"),"\t"))
            #else:
            #    foutPANCAN.write("\t"+string.join(string.split(line[:-1],"\t")[1:],"\t") )
        #foutPANCAN.write("\n")
                
        #data normalization per gene
        while 1:
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
                #if i==0:
                #    foutPANCAN.write(data[0])
                for i in range(1,len(data)):
                    if data[i]=="":
                        fout.write("\t")
                        #foutPANCAN.write("\t")
                        continue
                    fout.write("\t"+str(data[i]-average))
                    #foutPANCAN.write("\t"+str(data[i]-average))
                fout.write("\n")

            #foutPANCAN.write("\n")

    #keys.append("PANCAN")
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
        J['dataProducer']="UCSC Xena team"
        J["sample_type"]="tumor"
        J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
        J['domain']="TCGA"
        J['owner']="TCGA"
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
        J["groupTitle"]= "TCGA "+TCGAUtil.cancerGroupTitle[cancer]
        
        J['label']="gene expression RNAseq (IlluminaHiSeq pancan normalized)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq (IlluminaHiSeq), pancan normalized"
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq, mean-normalized (per gene) across all TCGA cohorts."
        J["description"]= J["description"] +"<br><br>For comparing data within this cohort, we recommend to use the \"gene expression RNAseq\" dataset. For comparing with other TCGA cohorts, we recommend to use the pancan normalized version of the \"gene expression RNAseq\" data. For comparing with data outside TCGA, we recommend using the percentile version if the non-TCGA data is normalized by percentile ranking."
        J["wrangling_procedure"]="Level_3 data (file names: *.rsem.genes.normalized_results) are download from TCGA DCC, log2(x+1) transformed, normalized across all TCGA cancer cohorts (all *.rsem.genes.normalized_results) per gene, and deposited into UCSC Xena repository."
        J['tags']=["cancer"] + TCGAUtil.tags[cancer]

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
                    origin.extend(value)

            J["tags"]= J["tags"]+["PANCAN12"]
            J["tags"]= J["tags"]+origin
            J["url"]="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
            J["name"]= "TCGA_PANCAN_exp_HiSeqV2_PANCAN"
        fout.write(json.dumps(J,indent=-1))
        fout.close()
    return
