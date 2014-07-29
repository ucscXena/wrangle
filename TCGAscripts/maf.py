import sys,string,os
import json,datetime
import math
import inspect
import copy

LEVEL="Level_2"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *
import xenaToMatrix
import radia

tmpDir="tmptmp/"

def unc_mixed_dnaseq_cont_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "Mixed_DNASeq_Cont_automated."
    PLATFORM = "mixed"
    suffix     = "unc"
    namesuffix = "mutation_unc_gene"
    dataProducer = "University of North Carolina"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def ucsc_illuminaga_dnaseq_cont (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_Cont."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    namesuffix = "mutation_ucsc_gene_maf"
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)
    
    clean=0
    namesuffix = "mutation_ucsc"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def ucsc_illuminaga_dnaseq_cont_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_Cont_automated."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_gene_maf"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean=0
    namesuffix = "mutation_ucsc"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)


def ucsc_solid_dnaseq_cont  (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "SOLiD_DNASeq_Cont."
    PLATFORM = "SOLiD"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_solid_gene_maf"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean=0
    namesuffix = "mutation_ucsc"
    radia.radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def broad_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def broad_illuminaga_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_automated"
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def broad_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "broad"
    namesuffix = "mutation_curated_broad_gene"
    dataProducer = "Broad Institute Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_curated_broad"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def hgsc_solid_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "SOLiD_DNASeq."
    PLATFORM = "SOLiD"
    suffix     = "bcm SOLiD"
    namesuffix = "mutation_bcm_solid_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_bcm_solid"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def hgsc_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)


def hgsc_illuminaga_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_automated"
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def hgsc_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "bcm"
    namesuffix = "mutation_curated_bcm_gene"
    dataProducer = "Baylor College of Medicine Human Genome Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_curated_bcm"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def wustl_illuminaga_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq."
    PLATFORM = "IlluminaGA"
    suffix     = "wustl"
    namesuffix = "mutation_wustl_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_wustl"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def wustl_illuminaga_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_automated"
    PLATFORM = "IlluminaGA"
    suffix     = "wustl"
    namesuffix = "mutation_wustl_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_wustl"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def wustl_illuminaga_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_curated"
    PLATFORM = "IlluminaGA"
    suffix     = "wustl"
    namesuffix = "mutation_curated_wustl_gene"
    dataProducer = "Genome Institute at Washington University Sequencing Center"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_curated_wustl"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def bcgsc_illuminahiseq_dnaseq (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_DNASeq."
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def bcgsc_illuminahiseq_dnaseq_automated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_DNASeq_automated"
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def bcgsc_illuminahiseq_dnaseq_curated (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_DNASeq_curated"
    PLATFORM = "IlluminaHiSeq"
    suffix     = "bcgsc"
    namesuffix = "mutation_curated_bcgsc_gene"
    dataProducer = "Michael Smith Genome Sciences Centre (British Columbia Genome Sciences Centre, BCGSC)"
    clean=1
    mafToMatrix (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

    clean =0
    namesuffix = "mutation_curated_bcgsc"
    mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def mafToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM):
    garbage=[tmpDir]
    os.system("rm -rf tmp_*") 
    if os.path.exists( tmpDir ):
        if clean:
            os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #multiple files in dir mode
    lastRelease={}
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

    #data processing multiple dirs mode
    if REALRUN:
        xena = outDir+cancer+"/"+cgFileName 
        fout= open(xena,'w')
        fout.write("#"+string.join(["sample","chr","start","end","gene","reference","alt","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =".maf"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    print infile
                    process_xena (infile, fout)
        fout.close()

    if not os.path.exists(outDir+cancer+"/"+cgFileName):
        return
    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")    
    J={}
    #stable    
    J["cgDataVersion"]=1
    J[":dataSubType"]="somatic mutation"
    J["redistribution"]= True
    J["dataProducer"]= dataProducer
    J["type"]= "mutationVector" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["PLATFORM"]= PLATFORM 
    if string.find( dataProducer ,"Broad")!=-1:
        J["method"]= "MutDect"
    elif string.find( dataProducer ,"Baylor")!=-1:
        J["method"]= "Baylor pipeline"
    elif string.find( dataProducer ,"Genome Institute at Washington University Sequencing Center")!=-1:
        J["method"]= "WashU pipeline"
    elif string.find( dataProducer ,"Michael Smith Genome Sciences Centre")!=-1:
        J["method"]= "BCGSC pipeline"
    elif string.find( dataProducer ,"University of Californis Santa Cruz GDAC")!=-1:
        J["method"]= "UCSC pipeline"
    elif string.find( dataProducer ,"University of North Carolina")!=-1:
        J["method"]= "UNC pipeline"
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .maf file from TCGA DCC, processed into UCSC Xena muation format into Xena repository"

    #change description
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["shortTitle"]= cancer +" mutation ("+suffix+" curated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+" curated)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["shortTitle"]= cancer +" mutation ("+suffix+" automated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+" automated)"
    else:
        J["shortTitle"]= cancer +" mutation ("+suffix+")"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+")"

    J["assembly"]="hg19"
    J["wholeGenome"]= True
    J["label"] = J["shortTitle"] 
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]
    J['owner']="TCGA"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data. Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using "+ J["method"] +" method."

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

def mafToMatrix (inDir, outDir, cancer,flog,PATHPATTERN,suffix, namesuffix, dataProducer,REALRUN, clean, PLATFORM):
    garbage=[tmpDir]
    os.system("rm -rf tmp_*") 
    if os.path.exists( tmpDir ):
        if clean:
            os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #multiple files in dir mode
    lastRelease={}
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

    #data processing multiple dirs mode
    if REALRUN:
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
                    process (infile, allGenes, samples, genes, dic)

        if len(samples)!=0:
            outputMatrix(samples, genes, dic, outDir+cancer+"/"+cgFileName)

    if not os.path.exists(outDir+cancer+"/"+cgFileName):
        return

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")    
    J={}
    #stable    
    J["cgDataVersion"]=1
    J[":dataSubType"]="somatic mutation"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J[":probeMap"]= "hugo"
    J["PLATFORM"]= PLATFORM 
    if string.find( dataProducer ,"Broad")!=-1:
        J["method"]= "MutDect"
    elif string.find( dataProducer ,"Baylor")!=-1:
        J["method"]= "Baylor pipeline"
    elif string.find( dataProducer ,"Genome Institute at Washington University Sequencing Center")!=-1:
        J["method"]= "WashU pipeline"
    elif string.find( dataProducer ,"Michael Smith Genome Sciences Centre")!=-1:
        J["method"]= "BCGSC pipeline"
    elif string.find( dataProducer ,"University of Californis Santa Cruz GDAC")!=-1:
        J["method"]= "UCSC pipeline"
    elif string.find( dataProducer ,"University of North Carolina")!=-1:
        J["method"]= "UNC pipeline"
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .maf file from TCGA DCC, processed into gene by sample matrix at UCSC into cgData repository"
    J["gain"]=10
    J["min"]=-0.1
    J["max"]=0.1

    #change
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+" curated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" curated)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+" automated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" automated)"
    else:
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+")"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+")"

    J["label"] = J["shortTitle"] 
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]
    J['owner']="TCGA"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data.  Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using "+ J["method"] +" method. <br><br> Red (=1) indicates that a non-silent somatic mutation (nonsense, missense, frame-shif indels, splice site mutations, stop codon readthroughs, change of start codon, inframe indels) was identified in the protein coding region of a gene, or any mutation identified in a non-coding gene. White (=0) indicates that none of the above mutation calls were made in this gene for the specific sample.<br><br>"
    J["description"] = J["description"] +"<br><br>"

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

def process (file, allGenes, samples, genes, dic):
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
                if header[i]=="t_alt_count":
                    t_alt_count =i
                if header[i]=="t_ref_count":
                    t_ref_count =i
                if header[i]=="Tumor_Sample_Barcode":
                    Tumor_Sample_Barcode =i

            continue

        #cut -f 1,9,16 brca_cleaned_filtered.maf |grep $GENE
        gene = data[Hugo_Symbol]
        try:
            mtype = xenaToMatrix.typeDic[data[Variant_Classification]]
        except:
            print "muation type not seen before", data[Variant_Classification]
            continue
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

def outputMatrix(samples, genes, dic, outfile):
    print len(genes), len(dic)
    fout=open(outfile,'w')
    fout.write("sample\t"+string.join(samples,"\t")+"\n")
    for gene in dic:
        fout.write(gene)
        for sample in samples:
            try:
                mtype = dic[gene][sample]
                if mtype < xenaToMatrix.non_silent_cutoff :
                    fout.write("\t1")
                else:
                    fout.write("\t0")
            except:
                fout.write("\t0")
        fout.write("\n")
    fout.close()

def process_xena (file, fout):
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
                if string.lower(header[i])=="t_ref_count":
                    t_ref_count =i
                if header[i]=="Tumor_Sample_Barcode":
                    Tumor_Sample_Barcode =i
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
            DNA_VAF = float(data[t_alt_count])/float(data[t_ref_count]+data[t_alt_count])
        except:
            DNA_VAF= ""
        RNA_VAF= ""
        AA_Change=""
        
        fout.write(string.join([sample,chrom,start,end, ref, alt, gene, mtype, str(DNA_VAF), str(RNA_VAF),AA_Change],"\t")+"\n")

    fin.close()
    
