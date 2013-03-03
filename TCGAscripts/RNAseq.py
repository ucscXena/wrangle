import sys,string,os
import json,datetime
import math
import inspect

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/illuminaga_rnaseq/rnaseq/
# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/illuminahiseq_rnaseq/rnaseq/
# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/bcgsc.ca/illuminaga_rnaseq/rnaseq/
# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/

def illuminahiseq_rnaseq_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_RNASeq"
    suffix = "IlluminaHiSeq_RNASeqV1"
    namesuffix = "HiSeqV1"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN)
    return

def illuminahiseq_rnaseqV2_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_RNASeqV2"
    suffix     = "IlluminaHiSeq_RNASeqV2"
    namesuffix = "HiSeqV2"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN)
    return

def illuminahiseq_rnaseqV2_exon_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_RNASeqV2"
    suffix     = "IlluminaHiSeq_RNASeqV2_exon"
    namesuffix = "HiSeqV2_exon"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN)
    return

def illuminaga_rnaseq_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_RNASeq"
    suffix     = "IlluminaGA_RNASeq"
    namesuffix = "GA"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN)
    return


def illuminaga_rnaseq_bcgsc (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN = "IlluminaGA_RNASeq"
    suffix      = "IlluminaGA_RNASeq"
    namesuffix = "GA"
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN)
    return

def geneRPKM (inDir, outDir, cancer,flog,PATHPATTERN,suffix, namesuffix, dataProducer,REALRUN):
    garbage=["tmptmp/"]
    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")

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
        if string.find(file,".tar.gz")!=-1 and REALRUN:
            os.system("tar -xzf "+inDir+file +" -C tmptmp/") 
            rootDir ="tmptmp/"
            
    #make sure there is data
    if REALRUN and (rootDir =="" or not os.path.exists(rootDir)):
        cleanGarbage(garbage)
        print "ERROR expect data, but wrong dirpath", rootDir, cancer, __name__
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= suffix

    #data processing multiple dirs mode
    if REALRUN:
        dataMatrix=[]
        samples={}
        genes={}

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample =""
                #v1
                pattern ="gene.quantification"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[1]
                    # bcgsc stupid sample name in file name
                    elif dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0]
                    else:
                        print "please check how to identify sample name"
                #v2
                pattern ="rsem.genes.normalized_results"
                if string.find(file,pattern)!=-1 and suffix== "IlluminaHiSeq_RNASeqV2":
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"
                #v2 exon from unc
                pattern ="bt.exon_quantification" 
                if string.find(file,pattern)!=-1 and suffix== "IlluminaHiSeq_RNASeqV2_exon":
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"

                if sample=="":
                    continue
                if sample in samples:
                    message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__
                    flog.write(message+"\n")
                    print message
                    continue
                # Test for barcode or UUID     #throw out all normals and control Analyte
                if sample[0:4]!="TCGA":
                    if TCGAUtil.UUID_CELLLINE.has_key(sample):
                        print "control cell line ignore", sample
                        continue
                else:
                    sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
                    if sampleTypeCode == False: # likely a uuid
                        continue
                    elif sampleTypeCode in ["20"]:
                        print "control cell line ignore", sample
                        continue

                p=len(samples)
                samples[sample]=p
            
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample=""
                #v1
                pattern ="gene.quantification"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[1]
                    # bcgsc stupid sample name in file name
                    elif dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0]
                    else:
                        print "please check how to identify sample name"
                    valuePOS=3
                    LOG2=1
                #v2
                pattern ="rsem.genes.normalized_results"
                if string.find(file,pattern)!=-1 and suffix== "IlluminaHiSeq_RNASeqV2":
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"
                    valuePOS=1
                    LOG2=1
                #v2 exon from unc
                pattern ="bt.exon_quantification" 
                if string.find(file,pattern)!=-1 and suffix== "IlluminaHiSeq_RNASeqV2_exon":
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"
                    valuePOS=3
                    LOG2=1
                    
                if sample=="":
                    continue
                if sample not in samples:
                    continue
                
                process(dataMatrix,samples,sample,genes, cancer,infile,flog, valuePOS,LOG2)
            
        outfile = outDir+cancer+"/"+cgFileName
        outputMatrix(dataMatrix, samples, genes, outfile, flog)
    
    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    J["cgDataVersion"]=1
    J[":dataSubType"]="geneExp"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["colNormalization"]=True
    J["PLATFORM"]= PATHPATTERN
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

    if suffix =="IlluminaHiSeq_RNASeqV2":  
        J["gain"]=0.67
    else:
        J["gain"]=1.0

    if PATHPATTERN in ["IlluminaHiSeq_RNASeq","IlluminaHiSeq_RNASeqV2"]:
        platformTitle ="Illumina HiSeq 2000 RNA Sequencing platform"
    if PATHPATTERN =="IlluminaGA_RNASeq":
        platformTitle =" Illumina Genome Analyzer RNA Sequencing platform"

    #change description
    J["description"]=""
    if string.find(suffix, "exon")==-1:
        J[":probeMap"]= "hugo"
        J["shortTitle"]="Gene Expression ("+suffix+")"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression ("+suffix+")"
        J["notes"]= "the probeMap should be tcgaGAF, but untill the probeMap is made, we will have to use hugo for the short term, however probably around 10% of the gene symbols are not HUGO names, but ENTRE genes"
        J["description"]= J["description"] +"The dataset shows TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression."+ \
                          " The gene expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
                          " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the gene-level transcription estimates, "
        
    else:
        J[":probeMap"]= "unc_RNAseq_exon"
        J["shortTitle"]="Exon Expression"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") exon expression ("+PATHPATTERN+")"
        J["description"]= J["description"] +"The dataset shows TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") exon expression."+ \
                          " The exon expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
                          " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the exon-level transcription estimates, "
        
    if suffix =="IlluminaHiSeq_RNASeqV2":
        J["description"] = J["description"] + "as in RSEM normalized count."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.rsem.genes.normalized_results) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"
    elif suffix =="IlluminaHiSeq_RNASeqV2_exon":
        J["description"] = J["description"] + "as in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.exon_quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"
    else:
        J["description"] = J["description"] + "as in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.gene.quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"

    if string.find(suffix, "exon")==-1:
        J["description"] = J["description"] + " Genes are mapped onto the human genome coordinates using UCSC cgData HUGO probeMap."
    else:
        J["description"] = J["description"] + " Exons are mapped onto the human genome coordinates using UCSC cgData unc_RNAseq_exon probeMap."
        
    if dataProducer =="University of North Carolina TCGA genome characterization center":
        J["description"] = J["description"] +\
                           " Reference to method description from "+dataProducer+": <a href=\"" + TCGAUtil.remoteBase +string.replace(inDir,TCGAUtil.localBase,"") +remoteDataDirExample+"/DESCRIPTION.txt\" target=\"_blank\"><u>DCC description</u></a>"
        
    J["description"] = J["description"] +\
                       "<br><br>In order to more easily view the differential expression between samples, we set the default view to center each gene or exon to zero by independently subtracting the mean of the genomic location on the fly. Users can view the original non-normalized values by uncheck the \"Normalize\" option. For more information on how to use the cancer browser, please refer to the help page."
    J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
    
    #change cgData
    J["name"]="TCGA_"+cancer+"_exp_"+namesuffix
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
            
    cleanGarbage(garbage)
    return

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf dir")
    return

def process(dataMatrix,samples, sample,genes, cancer,infile,flog, valuePOS, LOG2):
    # one sample a file
    fin=open(infile,'U')    
    fin.readline()

    l=[]
    for j in range (0,len(samples)):
        l.append("")
        
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        hugo = data[0]
        value= data[valuePOS]
        hugo = string.split(hugo,"|")[0]
        if hugo=="?":
            continue

        if hugo not in genes:
            p=len(genes)
            genes[hugo]=p
            dataMatrix.append(l[:])
        
        if value not in ["","null","NULL","Null","NA"]:
            if LOG2:
                value = float(value)
                if value<0:
                    value = ""
                else:
                    value = str(math.log10(float(value+1))/math.log10(2))
            x=genes[hugo]
            y=samples[sample]
            dataMatrix[x][y]=value
    fin.close()
    return 

def outputMatrix(dataMatrix, samples, genes, outfile, flog):
    fout = open(outfile,"w")
    fout.write("sample")
    for sample in samples:
        fout.write("\t"+sample)
    fout.write("\n")

    for gene in genes:
        fout.write(gene)
        for sample in samples:
            value = dataMatrix[genes[gene]][samples[sample]]
            if value !="":
                fout.write("\t"+value)
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
