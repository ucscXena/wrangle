import sys,string,os
import json,datetime
import math
import inspect
import copy
import uuid

LEVEL="Level_2"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *
import xenaToMatrix
import radia

tmpDir =str(uuid.uuid4())+"/"
VAF_cut =0.04

## REALRUN
# 2: only run matrix json
# 1: run both xena and matrix, real data rerun
# 0: run both json

#vcf protected germline
def germline_ucsc_illuminaga_dnaseq_cont (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution= False
    PATHPATTERN= "IlluminaGA_DNASeq_Cont."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "germline_ucsc_vcf"
    dataProducer = "University of Californis Santa Cruz GDAC"

    clean=0
    type ="germline"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution, type)

#vcf protected
def ucsc_illuminaga_dnaseq_cont_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution= False
    PATHPATTERN= "IlluminaGA_DNASeq_Cont_automated."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_maf_gene"
    dataProducer = "University of Californis Santa Cruz GDAC"

    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_ucsc_maf"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean=0
    namesuffix = "mutation_ucsc_vcf"
    type ="somatic"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution,type)

#maf open ucsc
def ucsc_illuminaga_dnaseq_automated (inDir, outDir, cancer, flog, REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq_automated."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_maf_gene"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean=0
    namesuffix = "mutation_ucsc_maf"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

#vcf
def ucsc_solid_dnaseq_cont  (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "SOLiD_DNASeq_Cont."
    PLATFORM = "SOLiD"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_solid_gene_maf"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean=0
    namesuffix = "mutation_ucsc"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def broad_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean=0
    namesuffix = "mutation_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def broad_illuminaga_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    VAF= True
    PATHPATTERN= "IlluminaGA_DNASeq_automated"
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution, VAF)

    clean =0
    namesuffix = "mutation_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution, VAF)

def broad_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_curated_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_curated_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def hgsc_solid_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "SOLiD_DNASeq_curated."
    PLATFORM = "SOLiD"
    suffix     = "bcm SOLiD"
    namesuffix = "mutation_curated_bcm_solid_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_curated_bcm_solid"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def hgsc_solid_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "SOLiD_DNASeq."
    PLATFORM = "SOLiD"
    suffix     = "bcm SOLiD"
    namesuffix = "mutation_bcm_solid_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_bcm_solid"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def hgsc_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)


def hgsc_illuminaga_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq_automated"
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def hgsc_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_curated_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_curated_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def wustl_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "wustl"
    namesuffix = "mutation_wustl_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_wustl"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def wustl_illuminahiseq_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaHiSeq_DNASeq_automated"
    PLATFORM = "IlluminaHiSeq"
    suffix     = "wustl_hiseq"
    namesuffix = "mutation_wustl_hiseq_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_wustl_hiseq"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def wustl_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "wustl"
    namesuffix = "mutation_curated_wustl_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_curated_wustl"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def bcgsc_illuminahiseq_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaHiSeq_DNASeq."
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def bcgsc_illuminahiseq_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaHiSeq_DNASeq_automated"
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def bcgsc_illuminahiseq_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    distribution = True
    PATHPATTERN= "IlluminaHiSeq_DNASeq_curated"
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_curated_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

    clean =0
    namesuffix = "mutation_curated_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM, distribution, VAF= False):
    #VAF =true is saying use VAF_cut, for broad data mostly 
    if REALRUN == 2:
        return

    REALRUN = 1 # always do real run 
    clean =1 #always clean up and rerun

    garbage=[tmpDir]
    os.system("rm -rf tmp_*")
    if os.path.exists( tmpDir ):
        if clean:
            os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #multiple files in dir mode
    lastRelease={}
    onlyArchive=0
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue

        if not os.path.exists(inDir +file+".md5"):
            print "file has no matching .md5 throw out", file
            continue

        #find lastest in each archive
        info = string.split(file,".")
        archive = info [-5]
        release = int(info [-4])

        if archive =="100":
            onlyArchive = "100"

        if not lastRelease.has_key(archive):
            lastRelease[archive]= release
        else:
            if lastRelease[archive]< release:
                lastRelease[archive]=release


    rootDir =""
    lastDate=None
    remoteDataDirExample =""

    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue

        if not os.path.exists(inDir +file+".md5"):
            continue

        #find the file that is the lastest release for the archive
        info = string.split(file,".")
        archive = info [-5]
        release = int(info [-4])

        if release != lastRelease[archive]:
            continue


        if onlyArchive and archive!=onlyArchive:
             continue


        #file latest date
        newDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)
        if not lastDate:
            lastDate = newDate
        if lastDate < newDate:
            lastDate = newDate

        if remoteDataDirExample =="":
            remoteDataDirExample = file[:-7]

        #is tar.gz?, uncompress multiple file mode
        if not clean:
            rootDir =tmpDir
        elif string.find(file,".tar.gz")!=-1 and REALRUN and clean:
            os.system("tar -xzf "+inDir+file +" -C "+tmpDir)
            rootDir =tmpDir

    #make sure there is data

    if REALRUN and (rootDir =="" or not os.path.exists(rootDir)):
        print "ERROR expect data, but wrong dirpath", rootDir, cancer, __name__
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= namesuffix
    assembly= None
    totalThrowout =0
    totalKept =0
    xena = outDir+cancer+"/"+cgFileName
    #data procesing multiple dirs mode
    if REALRUN:
        found =0
        fout= open(xena,'w')
        fout.write(string.join(["sample","chr","start","end","reference","alt","gene","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =".maf"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    print infile
                    found =1
                    data_assembly, throwout, kept = process_xena (infile, fout, assembly, VAF)
                    if not assembly:
                        assembly = data_assembly
                    if data_assembly != assembly:
                        continue 
                    totalKept = totalKept +kept
                    totalThrowout = totalThrowout +throwout
        fout.close()

        if not found:
            os.system("rm -f "+xena)

    if not os.path.exists(xena):
        os.system("rm -f "+xena+".json")
        return

    oHandle = open(xena+".json","w")
    J={}
    #stable
    J["dataSubType"]="somatic mutation (SNPs and small INDELs)"
    J["redistribution"]= distribution
    J["dataProducer"]= dataProducer
    J["type"]= "mutationVector"
    J["start_index"]=1
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["PLATFORM"]= PLATFORM
    if string.find( dataProducer ,"Broad")!=-1:
        J["method"]= "MuTect"
    elif string.find( dataProducer ,"Baylor")!=-1:
        J["method"]= "Baylor pipeline"
    elif string.find( dataProducer ,"Genome Institute at Washington University Sequencing Center")!=-1:
        J["method"]= "WashU pipeline"
    elif string.find( dataProducer ,"Michael Smith Genome Sciences Centre")!=-1:
        J["method"]= "BCGSC pipeline"
    elif string.find( dataProducer ,"University of Californis Santa Cruz GDAC")!=-1:
        J["method"]= "RADIA"
    elif string.find( dataProducer ,"University of North Carolina")!=-1:
        J["method"]= "UNC pipeline"
    else:
        J["method"]= ""

    #multiple dirs
    if string.find(inDir,TCGAUtil.localBase)!=-1:
        J["url"]=TCGAUtil.remoteBase \
            +string.replace(inDir,TCGAUtil.localBase,"")
    else:
        J["url"]=TCGAUtil.remoteBase \
            +string.replace(inDir,"/inside/depot/","")

    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .maf file from TCGA DCC, removed any calls from WGS samples, processed into UCSC Xena mutation format, stored in the UCSC Xena repository"

    #change description
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["label"]= suffix+" curated"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation ("+suffix+" curated)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["label"]= suffix+" automated"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation ("+suffix+" automated)"
    else:
        J["label"]= suffix
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation ("+suffix+")"

    J["assembly"]=assembly
    J["target"]= "whole exome"
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]
    J['owner']="TCGA"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data. Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using the "+ J["method"] +" method."
    if VAF:
        J["description"]= J["description"] + " When variant allele frequency information is available, only calls with VAF >"+ str(VAF_cut*100) +"% are kept, resulting in "+ str(totalKept) +" calls are kept and "+ str(totalThrowout) +" calls are removed."

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

def mafToMatrix (inDir, outDir, cancer,flog,PATHPATTERN,suffix, namesuffix, dataProducer,REALRUN, clean, PLATFORM, distribution, VAF=False):
    #VAF =true is saying use VAF_cut, for broad data mostly 
    garbage=[tmpDir]
    os.system("rm -rf tmp_*")
    if os.path.exists( tmpDir ):
        if clean:
            os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #multiple files in dir mode
    lastRelease={}
    onlyArchive=0

    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue

        if not os.path.exists(inDir +file+".md5"):
            print "file has no matching .md5 throw out", file
            continue

        #find lastest in each archive
        info = string.split(file,".")
        archive = info [-5]
        release = int(info [-4])

        if archive =="100":
            onlyArchive = "100"

        if not lastRelease.has_key(archive):
            lastRelease[archive]= release
        else:
            if lastRelease[archive]< release:
                lastRelease[archive]=release


    rootDir =""
    lastDate=None
    remoteDataDirExample =""
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue
        if not os.path.exists(inDir +file+".md5"):
            continue

        #find the file that is the lastest release for the archive
        info = string.split(file,".")
        archive = info [-5]
        release = int(info [-4])

        if release != lastRelease[archive]:
            continue

        if onlyArchive and archive!=onlyArchive:
            continue

        #file latest date
        newDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)
        if not lastDate:
            lastDate = newDate
        if lastDate < newDate:
            lastDate = newDate

        if remoteDataDirExample =="":
            remoteDataDirExample = file[:-7]

        #is tar.gz?, uncompress multiple file mode
        if not clean:
            rootDir =tmpDir
        elif string.find(file,".tar.gz")!=-1 and REALRUN and clean:
            os.system("tar -xzf "+inDir+file +" -C "+tmpDir)
            rootDir =tmpDir

    #make sure there is data
    if REALRUN ==1 and (rootDir =="" or not os.path.exists(rootDir)):
        print "ERROR expect data, but wrong dirpath", rootDir, cancer, __name__
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= namesuffix

    #data processing multiple dirs mode
    if REALRUN == 1:
        fgene=open("/data/TCGA/tcgaDataOneOff/Genomic/PANCAN/genes",'r')
        allGenes=[]
        for gene in string.split(fgene.read(),"\n"):
            allGenes.append(string.strip(gene))
        fgene.close()

        samples=[]
        genes=[]
        dic={}
        for gene in allGenes:
            dic[gene]={}

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =".maf"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    print infile
                    process (infile, allGenes, samples, genes, dic, VAF)

        if len(samples)!=0:
            xenaToMatrix.output(outDir+cancer+"/"+cgFileName, samples, genes, dic)

    if not os.path.exists(outDir+cancer+"/"+cgFileName):
        return

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    J={}
    #stable
    J["dataSubType"]="somatic non-silent mutation (gene-level)"
    J["redistribution"]= distribution
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J[":probeMap"]= "hugo"
    J["PLATFORM"]= PLATFORM
    J["unit"]= "binary non-silent mutation"
    if string.find( dataProducer ,"Broad")!=-1:
        J["method"]= "MuTect"
    elif string.find( dataProducer ,"Baylor")!=-1:
        J["method"]= "Baylor pipeline"
    elif string.find( dataProducer ,"Genome Institute at Washington University Sequencing Center")!=-1:
        J["method"]= "WashU pipeline"
    elif string.find( dataProducer ,"Michael Smith Genome Sciences Centre")!=-1:
        J["method"]= "BCGSC pipeline"
    elif string.find( dataProducer ,"University of Californis Santa Cruz GDAC")!=-1:
        J["method"]= "RADIA"
    elif string.find( dataProducer ,"University of North Carolina")!=-1:
        J["method"]= "UNC pipeline"
    else:
        J["method"]= ""

    #multiple dirs
    if string.find(inDir,TCGAUtil.localBase)!=-1:
        J["url"]=TCGAUtil.remoteBase \
            +string.replace(inDir,TCGAUtil.localBase,"")
    else:
        J["url"]=TCGAUtil.remoteBase \
            +string.replace(inDir,"/inside/depot/","")

    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "xena TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .maf file from TCGA DCC, processed into gene by sample matrix at UCSC, stored in the UCSC Xena repository"

    #change
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["label"]= suffix+" curated"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" curated)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["label"]= suffix+" automated"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" automated)"
    else:
        J["label"]= suffix
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+")"

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data.  Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using the "+ J["method"] +" method. <br><br> Red (=1) indicates that a non-silent somatic mutation (nonsense, missense, frame-shif indels, splice site mutations, stop codon readthroughs, change of start codon, inframe indels) was identified in the protein coding region of a gene, or any mutation identified in a non-coding gene. White (=0) indicates that none of the above mutation calls were made in this gene for the specific sample.<br><br>"

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

def process (file, allGenes, samples, genes, dic, VAF):
    fin =open(file, 'r')

    header=""
    for line in fin.readlines():
        if line[0] == "#":
            continue
        if string.strip(line)=="":
            continue
        data =string.split(line[:-1],'\t')
        if header=="":
            header=data
            for i in range (0,len(header)):
                if header[i]=="Hugo_Symbol":
                    Hugo_Symbol =i
                if header[i]=="NCBI_Build":
                    NCBI_Build =i
                if header[i]=="Variant_Classification":
                    Variant_Classification =i
                if header[i]=="Tumor_Sample_Barcode":
                    Tumor_Sample_Barcode =i
                if string.lower(header[i])=="t_alt_count":
                    t_alt_count =i
                if string.lower(header[i])=="t_ref_count":
                    t_ref_count =i
            continue

        gene = data[Hugo_Symbol]

        if gene =="":
            continue
        try:
            mtype = xenaToMatrix.typeDic[data[Variant_Classification]]
        except:
            print "mutation type not seen before", data[Variant_Classification]
            continue


        try:
            DNA_VAF = float(data[t_alt_count])/(float(data[t_ref_count])+float(data[t_alt_count]))
            if DNA_VAF >1.0:
                continue # crappy data

            if VAF: #require > VAF_cut to secure data accuracy
                if DNA_VAF <VAF_cut:
                    continue
        except:
            pass

        sample = data[Tumor_Sample_Barcode]

        if sample not in samples:
            samples.append(sample)
        if gene not in genes:
            genes.append(gene)
        if gene not in dic:
            dic[gene]={}
        if sample not in dic[gene]:
            dic[gene][sample]=mtype
        else:
            if mtype < dic[gene][sample]:
                dic[gene][sample]=mtype
    fin.close()

def process_xena (file, fout, assembly, VAF):
    fin =open(file, 'r')
    header=""
    ASSEMBLY=""
    throwout=0
    kept =0
    for line in fin.readlines():
        if line[0] == "#":
            continue
        if string.strip(line)=="":
            continue
        data =string.split(line[:-1],'\t')
        if header=="":
            header=data
            for i in range (0,len(header)):
                if header[i]=="Sequence_Source":
                    Sequence_Source =i
                if header[i]=="Hugo_Symbol":
                    Hugo_Symbol =i
                if header[i]=="NCBI_Build":
                    NCBI_Build =i
                if header[i]=="Chromosome":
                    Chromosome =i
                if string.lower(header[i])=="start_position":
                    Start_position =i
                if string.lower(header[i])=="end_position":
                    End_position =i
                if header[i]=="Variant_Classification":
                    Variant_Classification =i
                if header[i]=="Reference_Allele":
                    Reference_Allele =i
                if header[i]=="Tumor_Seq_Allele2":
                    Tumor_Seq_Allele2 =i
                if string.lower(header[i])=="t_alt_count":
                    t_alt_count =i
                if string.lower(header[i])=="tumor_alt_count":
                    t_alt_count =i
                if string.lower(header[i])=="t_ref_count":
                    t_ref_count =i
                if string.lower(header[i])=="tumor_ref_count":
                    t_ref_count =i
                if string.lower(header[i])=="rna_tumor_alt_count":
                    rna_t_alt_count =i
                if string.lower(header[i])=="rna_tumor_ref_count":
                    rna_t_ref_count =i
                if header[i]=="Tumor_Sample_Barcode":
                    Tumor_Sample_Barcode =i
                ## aa change
                if header[i]=="Protein_Change":
                    Protein_Change =i
                if header[i]=="AAChange":
                    Protein_Change =i
                if string.find(header[i],"amino_acid_change")!=-1:
                    Protein_Change =i
            continue

        #do not use  Sequence_Source=="WGS"
        if data[Sequence_Source]=="WGS":
            continue

        if ASSEMBLY =="":
            if data[NCBI_Build] in ["hg38","38","GRCh38","GRCh38-lite"]:
                ASSEMBLY = "hg38"
            if data[NCBI_Build] in ["hg19","37","GRCh37","GRCh37-lite"]:
                ASSEMBLY = "hg19"
            if data[NCBI_Build] in ["hg18","36","GRCh36","GRCh36-lite"]:
                ASSEMBLY = "hg18"
            print assembly, ASSEMBLY
            if assembly and assembly != ASSEMBLY:
                continue

        # use the maf input
        sample = data[Tumor_Sample_Barcode]
        chrom= data[Chromosome]
        if string.lower(chrom[0:2])!="ch":
            chrom= "chr"+chrom
        start= data[Start_position]
        end = data[End_position]
        ref =data[Reference_Allele]
        alt= data[Tumor_Seq_Allele2]
        gene = data[Hugo_Symbol]
        mtype = data[Variant_Classification]
        try:
            DNA_VAF = float(data[t_alt_count])/(float(data[t_ref_count])+float(data[t_alt_count]))
            if DNA_VAF >1.0:
                throwout =throwout+1
                continue #crappy data

            if VAF: #require > VAF_cut to secure data accuracy
                if DNA_VAF <VAF_cut:
                    throwout =throwout+1
                    continue
        except:
            DNA_VAF= ""
        try:
            RNA_VAF = float(data[rna_t_alt_count])/(float(data[rna_t_ref_count])+float(data[rna_t_alt_count]))
            if RNA_VAF >1.0:
                continue #crappy data
        except:
            RNA_VAF= ""
        try:
            AA_Change= data[Protein_Change]
        except:
            AA_Change= ""
        kept = kept + 1
        fout.write(string.join([sample,chrom,start,end, ref, alt, gene, mtype, str(DNA_VAF), str(RNA_VAF), AA_Change],"\t")+"\n")

    fin.close()
    return ASSEMBLY, throwout, kept
