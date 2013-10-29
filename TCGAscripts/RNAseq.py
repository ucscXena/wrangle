import sys,string,os
import json,datetime
import math
import inspect
import copy

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmptmp/"

# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/illuminaga_rnaseq/rnaseq/
# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/bcgsc.ca/illuminaga_rnaseqv2/rnaseqv2/
# /inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/

def illuminahiseq_rnaseqV2_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_RNASeqV2"
    suffix     = "IlluminaHiSeq"
    namesuffix = "HiSeqV2"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean = 1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)

    print cancer, "illuminahiseq_rnaseqV2_exon_unc"
    PATHPATTERN= "IlluminaHiSeq_RNASeqV2"
    suffix     = "IlluminaHiSeq"
    namesuffix = "HiSeqV2_exon"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean =0
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return

def illuminahiseq_rnaseqV2_exon_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaHiSeq_RNASeqV2"
    suffix     = "IlluminaHiSeq"
    namesuffix = "HiSeqV2_exon"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return


def illuminaga_rnaseqV2_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_RNASeqV2"
    suffix     = "IlluminaGA"
    namesuffix = "GAV2"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)

    print cancer, "illuminaga_rnaseqV2_exon_unc"
    PATHPATTERN= "IlluminaGA_RNASeqV2"
    suffix     = "IlluminaGA"
    namesuffix = "GAV2_exon"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean =0
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    
    return


def illuminaga_rnaseqV2_exon_unc (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_RNASeqV2"
    suffix     = "IlluminaGA"
    namesuffix = "GAV2_exon"
    dataProducer = "University of North Carolina TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return


def illuminaga_rnaseq_bcgsc (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN = "IlluminaGA_RNASeq"
    suffix      = "IlluminaGA"
    namesuffix = "GA"
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return

def illuminahiseq_rnaseq_bcgsc  (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN = "IlluminaHiSeq_RNASeq"
    suffix      = "IlluminaHiseq"
    namesuffix = "HiSeq"
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return



def geneRPKM (inDir, outDir, cancer,flog,PATHPATTERN,suffix, namesuffix, dataProducer,REALRUN,clean):
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
        allSamples={}

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample =""
                #v1
                pattern =".gene.quantification"
                #stupid hg18 and hg19 issues, ignore files with hg19 in its name
                if string.find(file,pattern)!=-1:
                    if string.find(file,"hg19")!=-1:
                        continue
                    infile = rootDir+dataDir+"/"+file
                    # bcgsc stupid sample name in file name
                    if dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0]
                    else:
                        print "please check how to identify sample name"
                #v2
                pattern ="rsem.genes.normalized_results"
                if string.find(file,pattern)!=-1  and string.find(namesuffix,"exon")==-1:
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"
                #v2 exon from unc
                pattern ="bt.exon_quantification" 
                if string.find(file,pattern)!=-1 and string.find(namesuffix,"exon")!=-1:
                    infile = rootDir+dataDir+"/"+file
                    # unc stupid sample name in file name
                    if dataProducer =="University of North Carolina TCGA genome characterization center":
                        sample = string.split(file,".")[2]
                    else:
                        print "please check how to identify sample name"

                if sample=="":
                    continue
                if sample in allSamples:
                    print len(allSamples)
                    message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__ +file
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

                p=len(allSamples)
                allSamples[sample]=p
                    
        c=0
        dataMatrix=[]
        tmpSamples={}
        genes={}
        oldgenes={}
        files=[]
        GOOD=1
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample=""
                #v1
                pattern =".gene.quantification"
                if string.find(file,pattern)!=-1:
                    if string.find(file,"hg19")!=-1:
                        continue
                    infile = rootDir+dataDir+"/"+file
                    # bcgsc stupid sample name in file name
                    if dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0]
                    else:
                        print "please check how to identify sample name"
                    valuePOS=3
                    LOG2=1
                #v2
                pattern ="rsem.genes.normalized_results"
                if string.find(file,pattern)!=-1 and string.find(namesuffix,"exon")==-1:
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
                if string.find(file,pattern)!=-1 and string.find(namesuffix,"exon")!=-1:
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
                if sample in tmpSamples: #duplicated samples
                    continue
                if sample not in allSamples:
                    continue

                p=len(tmpSamples)
                tmpSamples[sample]=p
                
                c=c+1
                #print c
                process(dataMatrix,tmpSamples,sample,genes, cancer,infile,flog, valuePOS,LOG2,250)
                if (c % 250)==0:
                    tmpout="tmp_"+ str(int(c/250.0))
                    r =outputMatrix(dataMatrix, tmpSamples, genes, oldgenes, tmpout, flog)
                    if r:
                        GOOD=0
                        
                    dataMatrix=[]
                    tmpSamples={}
                    oldgenes=copy.deepcopy(genes)
                    genes ={}
                    files.append(tmpout)
                    
        if (c % 250)!=0:
            tmpout= "tmp_"+ str(int(c/250.0)+1)
            files.append(tmpout)
            r= outputMatrix(dataMatrix, tmpSamples, genes, oldgenes,tmpout, flog)
            if r:
                GOOD=0
                
        #paste all together
        outfile = outDir+cancer+"/"+cgFileName
        if GOOD:
            os.system("paste -d \'\' "+string.join(files," ")+" > "+ outfile)
        for file in files: 
            os.system("rm "+ file) 
        if not GOOD:
            sys.exit()
            
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
    J["gain"]=1.0

    if PATHPATTERN in ["IlluminaHiSeq_RNASeq","IlluminaHiSeq_RNASeqV2"]:
        platformTitle ="Illumina HiSeq 2000 RNA Sequencing platform"
    if PATHPATTERN in [ "IlluminaGA_RNASeq", "IlluminaGA_RNASeqV2"]:
        platformTitle =" Illumina Genome Analyzer RNA Sequencing platform"

    #change description
    J["description"]=""
    if string.find(namesuffix, "exon")==-1:
        J[":probeMap"]= "hugo"

        if cancer != "OV":
            J["shortTitle"]= cancer +" "+"gene expression ("+suffix+")"
            J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq ("+suffix+")"
        else:
            if dataProducer =="University of North Carolina TCGA genome characterization center":
                J["shortTitle"]= cancer +" "+"gene expression ("+suffix+" UNC)"
                J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq ("+suffix+" UNC)"
            else:
                J["shortTitle"]= cancer +" "+"gene expression ("+suffix+" BC)"
                J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq ("+suffix+" BC)"
                
        J["notes"]= "the probeMap should be tcgaGAF, but untill the probeMap is made, we will have to use hugo for the short term, however probably around 10% of the gene symbols are not HUGO names, but ENTRE genes"
        J["description"]= J["description"] +"TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression by RNAseq.<br><br>"+ \
                          " The gene expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
                          " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the gene-level transcription estimates, "
    else:
        J[":probeMap"]= "unc_RNAseq_exon"

        J["shortTitle"]= cancer +" "+"exon expression ("+suffix+")"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") exon expression by RNAseq ("+suffix+")"
        J["description"]= J["description"] +"TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") exon expression by RNAseq.<br><br>"+ \
                          " The exon expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
                          " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the exon-level transcription estimates, "
        
    if PATHPATTERN in [ "IlluminaHiSeq_RNASeqV2","IlluminaGA_RNASeqV2"] and string.find(namesuffix, "exon")==-1:
        J["description"] = J["description"] + "as in RSEM normalized count."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.rsem.genes.normalized_results) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"
    elif PATHPATTERN in [ "IlluminaHiSeq_RNASeqV2","IlluminaGA_RNASeqV2"] and string.find(namesuffix, "exon")!=-1:
        J["description"] = J["description"] + "as in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.exon_quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"
    else:
        J["description"] = J["description"] + "as in RPKM values (Reads Per Kilobase of exon model per Million mapped reads)."
        J["wrangling_procedure"]= "Level_3 Data (file names: *.gene.quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into cgData repository"

    if string.find(namesuffix, "exon")==-1:
        J["description"] = J["description"] + " Genes are mapped onto the human genome coordinates using UCSC cgData HUGO probeMap."
    else:
        J["description"] = J["description"] + " Exons are mapped onto the human genome coordinates using UCSC cgData unc_RNAseq_exon probeMap."
        
    if dataProducer =="University of North Carolina TCGA genome characterization center":
        J["description"] = J["description"] +\
                           " Reference to method description from "+dataProducer+": <a href=\"" + TCGAUtil.remoteBase +string.replace(inDir,TCGAUtil.localBase,"") +remoteDataDirExample+"/DESCRIPTION.txt\" target=\"_blank\"><u>DCC description</u></a>"
        
    J["description"] = J["description"] +\
                       "<br><br>In order to more easily view the differential expression between samples, we set the default view to center each gene or exon to zero by independently subtracting the mean of the genomic location on the fly. Users can view the original non-normalized values by uncheck the \"Normalize\" option. For more information on how to use the cancer browser, please refer to the help page."
    #J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc

    J["label"] = J["shortTitle"] 
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]
    J['owner']="TCGA"
    J['gdata_tags'] =["transcription"]
    
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
    
    return

def process(dataMatrix,samples, sample,genes, cancer,infile,flog, valuePOS, LOG2, maxLength):
    # one sample a file
    fin=open(infile,'U')    
    fin.readline()

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
            l=[]
            for j in range (0,maxLength):
                l.append("")    
            dataMatrix.append(l)
        
        if value not in ["","null","NULL","Null","NA"]:
            if LOG2:
                value = float(value)
                if value<0:
                    value = ""
                else:
                    value = math.log10(float(value+1))/math.log10(2)

            x=genes[hugo]
            y=samples[sample]
            dataMatrix[x][y]=value
            
    fin.close()
    return 

def outputMatrix(dataMatrix, samples, genes, oldgenes,outfile, flog):
    #compare genes and oldgenes:
    if oldgenes!={}:
        if len(genes)!=len(oldgenes):
            print "ERROR genes total length is different"
            return 1
        for gene in genes:
            if genes[gene]!=oldgenes[gene]:
                print "ERROR gene order is different",gene
                return 1
            
    fout = open(outfile,"w")
    if oldgenes=={}:
        fout.write("sample")
    for sample in samples:
        fout.write("\t"+sample)
    fout.write("\n")

    for gene in genes:
        if oldgenes=={}:
            fout.write(gene)
        for sample in samples:
            value = dataMatrix[genes[gene]][samples[sample]]
            if value !="":
                value = "%.4f" % (float(value))
                fout.write("\t"+value)
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
    return 0
