import string, os, sys
import json,datetime
import csv

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *
from CGDataUtil import *
from CGDataLib import *
import TCGAUtil
import xenaToMatrix

#for individual
def RNAseq (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return    
    print cancer, sys._getframe().f_code.co_name
    filename = "HiSeqV2"
    processRNA(filename, dir,outDir, cancer,flog, REALRUN)

#for combination
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
    memVersion = True
    processMatrix (filename, dir, outDir, cancer,flog, REALRUN, memVersion)

def Methylation27 ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "HumanMethylation27"
    memVersion = True
    processMatrix (filename, dir, outDir, cancer,flog, REALRUN, memVersion)

def Methylation450 ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "HumanMethylation450"
    memVersion = True
    processMatrix (filename, dir, outDir, cancer,flog, REALRUN, memVersion)

def miRNA  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name

    filename = "miRNA_HiSeq"
    processMatrix (filename, dir, outDir, cancer,flog, REALRUN)

    filename = "miRNA_HiSeq.probeMap"
    processProbeMap (filename, dir, outDir, cancer, REALRUN)

    filename = "miRNA_GA"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

    filename = "miRNA_GA.probeMap"
    processProbeMap (filename, dir, outDir, cancer, REALRUN)

    filename = "miRNA_HiSeq_gene"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

    filename = "miRNA_GA_gene"
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

def mutation_ucsc_maf  ( dir, outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "mutation_ucsc_maf"
    processMutationData (filename, dir,outDir, cancer,flog, REALRUN)

def Gistic2 (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "Gistic2_CopyNumber_Gistic2_all_data_by_genes"
    processMatrix (filename, dir,outDir, cancer,flog, REALRUN)

def Gistic2_threshold (dir,outDir, cancer,flog,REALRUN):
    if cancer !="PANCAN":
        return
    print cancer, sys._getframe().f_code.co_name
    filename = "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
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

    features=["sample_type","sample_type_id","gender","age_at_initial_pathologic_diagnosis",
              "_OS","_OS_IND","_OS_UNIT",
              "_RFS","_RFS_IND","_RFS_UNIT",
              "_EVENT","_TIME_TO_EVENT","_TIME_TO_EVENT_UNIT",
              "_primary_site","_primary_disease"]   

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
                cfout.write("feature\tattribute\tvalue\n")
                featureConfig=1

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
        if cancer in ["FPPP", "LUNG","COADREAD","PANCAN","GBMLGG","PANCAN12"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename

        #if not os.path.exists(cancerFile):
        #    if filename =="miRNA":
        #        cancerFile = cancerDir+"/"+filename+"_HiSeq"
        #        if not os.path.exists(cancerFile):
        #            cancerFile = cancerDir+"/"+filename+"_GA"
        #            if not os.path.exists(cancerFile):
        #                continue
        #            else:
        #                inFiles.append(cancerFile)
        #        else:
        #            inFiles.append(cancerFile)
        #    else:
        #        continue

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
    elif filename=="mutation_ucsc_maf":
        mutation_ucsc_mafJSON(J,cancer,genelevel)
    elif filename=="mutation_bcgsc":
        mutation_bcgscJSON(J,cancer,genelevel)
    elif filename=="mutation_wustl":
        mutation_wustlSON(J,cancer,genelevel)
    commonJSON(J, cancer)
    J["dataSubType"]="somatic mutation (SNPs and small INDELs)"
    J["assembly"]="hg19"
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
    elif filename=="mutation_ucsc_maf":
        mutation_ucsc_mafJSON(J,cancer,genelevel)
    elif filename=="mutation_bcgsc":
        mutation_bcgscJSON(J,cancer,genelevel)
    elif filename=="mutation_wustl":
        mutation_wustlSON(J,cancer,genelevel)
    J[":probeMap"]="hugo"
    J["unit"] = "binary non-silent mutation"
    J["dataSubType"]="somatic non-silent mutation (gene-level)"
    commonJSON(J, cancer)

    fout.write(json.dumps(J,indent=-1))
    fout.close()

    return

def processProbeMap (filename, dir, outDir, cancer, REALRUN):
    inFiles =processFiles (filename, dir, outDir, cancer )

    cancer="PANCAN"
    outfile = outDir +cancer + "/"+filename

    if REALRUN:
        root = "/data/TCGA/"
        os.system("python mergeTableFiles.py.py noDup "+outfile+" "+root+" " + string.join(inFiles,' '))
    
    J={}
    J["type"]="probeMap"
    fout = open(outfile+".json","w")
    fout.close()
    return

def processMatrix (filename, dir,outDir, cancer,flog, REALRUN, memVersion = False):
    inFiles =processFiles (filename, dir,outDir, cancer )

    cancer="PANCAN"
    outfile = outDir +cancer + "/"+filename

    if REALRUN:
        root = "/data/TCGA/"
        os.system("python mergeGenomicMatrixFiles_memEfficient.py "+outfile+" "+root+" " + string.join(inFiles,' '))

    J={}
    fout = open(outfile+".json","w")

    cancer="PANCAN"    
    J["type"]= "genomicMatrix"
    J['url']="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
    if filename=="Gistic2_CopyNumber_Gistic2_all_data_by_genes":
        gisticJSON(J,cancer)
    elif filename=="Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes":
        gistic_thresholdJSON(J,cancer)
    elif filename=="HiSeqV2":
        HiSeqV2JSON(J,cancer)
    elif filename=="HiSeqV2_exon":
        HiSeqV2_exonJSON(J,cancer)
    elif filename=="miRNA_HiSeq":
        miRNA_HiSeqJSON(J,cancer)
    elif filename=="miRNA_GA":
        miRNA_GAJSON(J,cancer)
    elif filename=="miRNA_HiSeq_gene":
        miRNA_HiSeq_geneJSON(J,cancer)
    elif filename=="miRNA_GA_gene":
        miRNA_GA_geneJSON(J,cancer)
    elif filename=="HumanMethylation27":
        HumanMethylation27JSON(J,cancer)
    elif filename=="HumanMethylation450":
        HumanMethylation450JSON(J,cancer)
 
    commonJSON(J, cancer)

    fout.write(json.dumps(J,indent=-1))
    fout.close()
    return

def commonJSON(J, cancer):
    J['dataProducer']="UCSC Xena team"
    J["sample_type"]=["tumor"]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['tags']=["cancer"]+TCGAUtil.tags[cancer]
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
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
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation dataa, compiled using data from all TCGA cohorts all TCGA cohorts. The calls are generated at Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC) using the BCGSC pipeline method. BCGSC's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    

def mutation_wustlSON (J, cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_wustl"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (wustl)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" somatic mutation (wustl)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (wustl)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+genelevel+" nonsilent somatic mutation (wustl)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation dataa, compiled using data from all TCGA cohorts. The calls are generated at Washington University Genome Center using the WashU pipeline method. BCGSC's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"

def mutation_ucsc_mafJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_ucsc_maf"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (ucsc)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" somatic mutation (ucsc)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (ucsc)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" nonsilent somatic mutation (ucsc)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation dataa, compiled using data from all TCGA cohorts.  The calls are generated at University of Californis Santa Cruz GDAC using the RADIA method. RADIA's calls (mafs) from various TCGA cohorts are combined to produce this dataset. Reference to RADIA: PMID: 25405470."
    J["description"] = J["description"] +"<br><br>"    

def mutation_bcmJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_bcm"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (bcm)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+ " somatic mutation (bcm)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (bcm)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+ " nonsilent somatic mutation (bcm)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation dataa, compiled using data from all TCGA cohorts.  The calls are generated at Baylor College of Medicine Human Genome Sequencing Center using the Baylor pipeline method. Baylor's calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    

def mutation_broadJSON(J,cancer,genelevel):
    J['name']= "TCGA_PANCAN_mutation_broad"+string.strip(genelevel).replace("-","")
    if genelevel =="":
        J["label"]= "somatic mutation SNPs and small INDELs (broad)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" somatic mutation (broad)"
    else:
        J["label"]= "somatic gene-level non-silent mutation (broad)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+ genelevel+" nonsilent somatic mutation (broad)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" somatic mutation data, compiled using data from all TCGA cohorts. The calls are generated at Broad Institute Genome Sequencing Center using the MuTect method. MuTect calls from various TCGA cohorts are combined to produce this dataset."
    J["description"] = J["description"] +"<br><br>"    
    
def miRNA_HiSeqJSON(J,cancer):
    J['name']= "TCGA_PANCAN_miRNA_HiSeq"
    J["label"]= "miRNA isoform expression (HiSeq)"
    J["probeMap"] = "miRNA_HiSeq.probeMap"
    J["dataSubType"]="miRNA isoform expression RNAseq"
    J["idSubType"] = "isoform" 
    J["PLATFORM"] = "IlluminaHiSeq_miRNASeq"
    J["colNormalization"]=True
    J["unit"]="log2(RPM+1)"
    return

def miRNA_GAJSON(J,cancer):
    J['name']= "TCGA_PANCAN_miRNA_GA"
    J["label"]= "miRNA isoform expression (GA)"
    J["probeMap"] = "miRNA_GA.probeMap"
    J["dataSubType"]="miRNA isoform expression RNAseq"
    J["idSubType"] = "isoform"
    J["PLATFORM"] = "IlluminaGA_miRNASeq"
    J["colNormalization"]=True
    J["unit"]="log2(RPM+1)"
    return

def miRNA_HiSeq_geneJSON(J,cancer):
    J['name']= "TCGA_PANCAN_miRNA_HiSeq_gene"
    J["label"]= "miRNA gene expression (HiSeq)"
    if J.has_key(":probeMap"):
        J.pop(":probeMap")
    J["dataSubType"]="miRNA gene expression RNAseq"
    J["idSubType"] = "gene"
    J["PLATFORM"] = "IlluminaHiSeq_miRNASeq"
    J["colNormalization"]=True
    J["unit"]="log2(RPM+1)"
    return

def miRNA_GA_geneJSON(J,cancer):
    J['name']= "TCGA_PANCAN_miRNA_GA_gene"
    J["label"]= "miRNA gene expression (GA)"
    if J.has_key(":probeMap"):
        J.pop(":probeMap")
    J["dataSubType"]="miRNA gene expression RNAseq"
    J["idSubType"] = "gene"
    J["PLATFORM"] = "IlluminaGA_miRNASeq"
    J["colNormalization"]=True
    J["unit"]="log2(RPM+1)"
    return

def HumanMethylation27JSON(J,cancer):
    J['name']= "TCGA_PANCAN_HumanMethylation27"
    J["label"]= "DNA methylation (Methylation27K)"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" DNA methylation (HumanMethylation27K)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" DNA methylation 27K array beta values, compiled by combining available data from all TCGA cohorts. DNA methylation profile was measured experimentally using the Illumina Infinium HumanMethylation27 platform."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="illuminaMethyl27K_hg18_gpl8490"
    J["dataSubType"]="DNA methylation"
    J["unit"]="beta value"
    return

def HumanMethylation450JSON(J,cancer):
    J['name']= "TCGA_PANCAN_HumanMethylation450"
    J["label"]= "DNA methylation (Methylation450K)"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" DNA methylation (HumanMethylation450K)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" DNA methylation 450K array beta values, compiled by combining available data from all TCGA cohorts. DNA methylation profile was measured experimentally using the Illumina Infinium HumanMethylation450 platform."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="illuminaMethyl450_hg19_GPL16304"
    J["dataSubType"]="DNA methylation"
    J["unit"]="beta value"
    return

def HiSeqV2JSON (J, cancer):
    J['name']= "HiSeqV2_PANCAN"
    J["label"]= "gene expression"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" gene expression (polyA+ IlluminaHiSeq)"
    J["RNAtype"]="polyA+",
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" gene expression by RNAseq, compiled using data from all TCGA cohorts. Gene expression was measured using the IlluminaHiSeq technology. Data from all TCGA cohorts are combined to produce this dataset. Values are log2(x+1) transformed RSEM values.<br>"
    J[":probeMap"]="hugo"
    J["dataSubType"]="gene expression RNAseq"
    J["colNormalization"]=True
    J["unit"]="log2(normalized_count+1)"
    J["wrangling_procedure"] = "Level_3 data (file names: *.rsem.genes.normalized_results) are downloaded from each cancer project at TCGA DCC, log2(normalized_count+1) transformed, and then combined at UCSC into Xena repository."
    return

def HiSeqV2_exonJSON (J, cancer):
    J['name']= "HiSeqV2_exon_PANCAN"
    J["shortTitle"]= "exon expression"
    J["label"] = J["shortTitle"] 
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" exon expression (polyA+ IlluminaHiSeq)"
    J["RNAtype"]="polyA+",
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" exon expression by RNAseq, compiled using data from all TCGA cohorts. Exon expression was measured using the IlluminaHiSeq technology. Data from all TCGA cohorts are combined to produce this dataset. Values are log2(RPKM+1) transformed exon-level transcription estimates in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="unc_RNAseq_exon.hg19"
    J["dataSubType"]="exon expression RNAseq"
    J["colNormalization"]=True
    J["unit"]="log2(RPKM+1)"
    J["wrangling_procedure"] = "Level_3 data (file names: *.exon_quantification.txt) are downloaded from each cancer project at TCGA DCC, log2(x+1) transformed, and then combined at UCSC into Xena repository."
    return

def gisticJSON(J,cancer):
    J['name']= "TCGA_PANCAN_gistic2"
    J["label"]= "gene-level copy number (gistic2)"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" gene-level copy number (gistic2)"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" gene-level copy number variation (CNV) estimated using the GISTIC2 method, compiled using data from all TCGA cohorts. Copy number was measured experimentally using whole genome microarray at Broad TCGA genome characterization center. Subsequently, TCGA FIREHOSE pipeline applied GISTIC2 method to produce segmented CNV data, which was then mapped to genes to produce gene-level estimates. Gistic2 data from all TCGA cohorts are combined to produce this dataset. Reference to GISTIC2 method PMID:21527027."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="hugo"
    J["unit"] = "Gistic2 copy number"
    J["dataSubType"]="copy number (gene-level)"
    return

def gistic_thresholdJSON(J,cancer):
    J['name']= "TCGA_PANCAN_gistic2_threshold"
    J["label"]= "gene-level copy number (gistic2_thresholded)"
    J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" gene-level copy number (gistic2_thresholded)" 
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" gene-level copy number variation (CNV) estimated using the GISTIC2 threshold method, compiled using data from all TCGA cohorts. Copy number was measured experimentally using whole genome microarray at a TCGA genome characterization center. Subsequently, GISTIC2 method was applied using the TCGA FIREHOSE pipeline to produce gene-level copy number estimates. GISTIC2 further thresholded the estimated values to -2,-1,0,1,2, representing homozygous deletion, single copy deletion, diploid normal copy, low-level copy number amplification, or high-level copy number amplification. Genes are mapped onto the human genome coordinates using UCSC cgData HUGO probeMap. Reference to GISTIC2 method PMID:21527027."
    J["description"] = J["description"] +"<br><br>"    
    J[":probeMap"]="hugo"
    J["unit"] = "-2,-1,0,1,2: 2 copy del, 1 copy del, no change, amplification, high-amplification"
    J["dataSubType"]="copy number (gene-level)"
    return

def processRNA (filename, dir,outDir, cancer,flog, REALRUN):
    inFiles ={}
    outFiles={}
    for cancer in os.listdir(outDir):
        if cancer in ["PANCAN"]:
            continue

        cancerDir= outDir+ cancer
        cancerFile = cancerDir+"/"+filename
        if not os.path.exists(cancerFile):
            continue

        cancerOutFile = outDir+"/"+cancer+"/"+filename+"_PANCAN"

        inFiles[cancer]= cancerFile
        outFiles[cancer]= cancerOutFile

    keys = inFiles.keys()

    if REALRUN:
        #header:
        for i in range (0,len(keys)):
            fin = open(inFiles[keys[i]],'r')
            fout= open(outFiles[keys[i]],'w')
            inFiles[keys[i]] =fin
            outFiles[keys[i]]=fout
            line = fin.readline()
            fout.write(line)

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
                for i in range(1,len(data)):
                    if data[i]=="":
                        fout.write("\t")
                        continue
                    fout.write("\t"+str(data[i]-average))
                fout.write("\n")

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
        
        J["RNAtype"]="ployA+"
        J['name']= J['name']+"_PANCAN"
        J['dataProducer']="UCSC Xena team"
        J["sample_type"]="tumor"
        J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
        J['label']="gene expression RNAseq ("+ J["RNAtype"]+ " IlluminaHiSeq pancan normalized)"
        J['longTitle']="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq ("+ J["RNAtype"]+" IlluminaHiSeq), pancan normalized"
        J["unit"]="pan-cancer normalized log2(norm_counts+1)"
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq, mean-normalized (per gene) across all TCGA cohorts. Values in this dataset are generated at UCSC by combining \"gene expression RNAseq\" values of all TCGA cohorts, values are then mean-centered per gene, then extracting the converted data only belongs to the this cohort."

        J["description"]= J["description"] +"<br><br>For comparing data within this cohort, we recommend to use the \"gene expression RNAseq\" dataset. For questions regarding the gene expression of this particular cohort in relation to other types tumors, you can use the pancan normalized version of the \"gene expression RNAseq\" data. For comparing with data outside TCGA, we recommend using the percentile version if the non-TCGA data is normalized by percentile ranking. For more information, please see our Data FAQ: <a href=https://docs.google.com/document/d/1q-7Tkzd7pci4Rz-_IswASRMRzYrbgx1FTTfAWOyHbmk/edit?usp=sharing target=\"_blank\"><u>here</u></a>.<br>"

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
