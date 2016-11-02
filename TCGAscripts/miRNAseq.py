import sys,string,os
import json,datetime
import math
import inspect
import copy
import Jing_util
import mergeGenomicMatrixFiles_memEfficient as mergeGenomicMatrixFiles

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *
from processRSEM2percentile import *

tmpDir="tmpmiRNA/"

# /data/TCGA/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/bcgsc.ca/illuminaga_mirnaseq/mirnaseq/
# /data/TCGA/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/uvm/cgcc/bcgsc.ca/illuminahiseq_mirnaseq/mirnaseq/

def illuminaga_miRnaseq_bcgsc (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN = "IlluminaGA_miRNASeq"
    suffix      = "IlluminaGA"
    namesuffix = "miRNA_GA"
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer, REALRUN, clean)
    return

def illuminahiseq_miRnaseq_bcgsc  (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN = "IlluminaHiSeq_miRNASeq"
    suffix      = "IlluminaHiseq"
    namesuffix = "miRNA_HiSeq"
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"
    clean =1
    geneRPKM (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN, clean)
    return

def mergeGA_Hiseq (inDir, outDir, cancer, flog, REALRUN):
    print cancer, sys._getframe().f_code.co_name
    dataProducer = "British Columbia Cancer Agency TCGA genome characterization center"

    if not os.path.exists( outDir ):
        return
    file1 = outDir+cancer+"/"+"miRNA_GA"
    file2 = outDir+cancer+"/"+"miRNA_HiSeq"
    if ( (not os.path.exists(file1 )) or (not os.path.exists(file2 ))):
        return

    outfile = outDir+cancer+"/"+"miRNA"
    if REALRUN:
        genes={}
        samples={}
        dataMatrix=[]
        mergeGenomicMatrixFiles.header (samples, file1)
        mergeGenomicMatrixFiles.header (samples, file2)
        mergeGenomicMatrixFiles.process(genes, samples, dataMatrix, file1)
        mergeGenomicMatrixFiles.process(genes, samples, dataMatrix, file2)
        mergeGenomicMatrixFiles.outputMatrix(dataMatrix, samples, genes, outfile)

    #json
    J={}

    J["cgDataVersion"]=1
    J["redistribution"]= True
    J["dataProducer"]= dataProducer
    J["colNormalization"]=True
    J["PLATFORM"]= "IlluminaGA and IlluminaHiSeq miRNAseq"
    J["type"]= "genomicMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["unit"]="log2(RPM+1)"

    platformTitle ="Illumina Genome Analyzer and HiSeq 2000 RNA Sequencing platform"

    #change description
    J["description"]=""
    J["dataSubType"]="miRNA expression RNAseq"
    J["label"]= "miRNA isoform expression (GA, HiSeq)"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") miRNA expression by RNAseq (Illumina GA, HiSeq)"
    J["description"]= J["description"] +"TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") miRNA expression by RNAseq.<br><br>"+ \
        " The miRNA expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
        " Level 3 interpreted level data was downloaded from TCGA data coordination center. Download data is in unit reads per million mapped reads (RPM). This dataset shows the miRNA transcription estimates in log2 (RPM). For more information, see : http://nar.oxfordjournals.org/content/early/2015/08/13/nar.gkv808.full ."

    J["description"] = J["description"] +\
                       "<br><br>In order to more easily view the differential expression between samples, we set the default view to center each miRNA to zero by independently subtracting the mean of each miRNA across the cohort on the fly. Users can view the original non-normalized values by adjusting visualization settings."
    J["description"] = J["description"] +"<br><br>"

    J["wrangling_procedure"]= "Level_3 Data (file names: *.hg19.mirna.quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into Xena repository."

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    J['gdata_tags'] =["transcription","miRNA"]
    J["name"]="TCGA_"+cancer + "_miRNA"
    oHandle = open(outfile+".json",'w')
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()

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
        #hg19 or not
        pattern =".hg19.mirbase20.isoform.quantification"
        HG19= 0
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                if string.find(file,pattern)!=-1 :
                    HG19=1
                    break
            if HG19:
                break
        if HG19:
            pattern =".hg19.mirbase20.isoform.quantification"
        else:
            pattern =".isoform.quantification"

        allSamples={}
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample =""
                #bcgsc gene
                if string.find(file,pattern)!=-1 :
                    infile = rootDir+dataDir+"/"+file
                    # bcgsc stupid sample name in file name
                    if dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0][:15]
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
        mapping = {}
        tmpSamples={}
        genes={}
        oldgenes={}
        files=[]
        GOOD=1
        BATCH = 2500
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample=""
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    # bcgsc stupid sample name in file name
                    if dataProducer=="British Columbia Cancer Agency TCGA genome characterization center":
                        sample = string.split(file,".")[0][:15]
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
                process(dataMatrix, mapping, tmpSamples,sample, genes, cancer, infile,flog, valuePOS,LOG2,BATCH)

                if (c % BATCH)==0:
                    tmpout="tmp_"+ str(int(c/float(BATCH)))
                    r =outputMatrix(dataMatrix, tmpSamples, genes, oldgenes, tmpout, flog)
                    if r:
                        GOOD=0

                    dataMatrix=[]
                    tmpSamples={}
                    oldgenes=copy.deepcopy(genes)
                    genes ={}
                    files.append(tmpout)

        if (c % BATCH)!=0:
            tmpout= "tmp_"+ str(int(c/float(BATCH))+1)
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

        #probeMap
        probefile= outDir+cancer+"/"+namesuffix+".probeMap"
        outputProbeMap (probefile, mapping)


        #transcript datafile
        datafile= outDir+cancer+"/"+cgFileName
        if not os.path.exists(datafile):
            return

        #gene datafile
        genefile = outDir+cancer+"/"+cgFileName+"_gene"
        os.system("python ../support/genomicMatrixToGeneMatrix_memInEfficient.py "+ datafile + " " + probefile +' ' +  genefile + " add log2 1")

    #probeMap json
    probefile= outDir+cancer+"/"+namesuffix+".probeMap"
    fout = open(probefile +".json","w")
    J = {}
    J['type'] = 'probeMap'
    J['assembly'] ='hg19'
    fout.write( json.dumps( J, indent=-1 ) )
    fout.close()

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    J={}
    #stable
    J["cgDataVersion"]=1
    J["redistribution"]= True
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
    J["unit"]="log2(RPM+1)"
    J["idSubType"] = "isoform"
    J[":probeMap"]= "miRNA_isoform_TCGA.probeMap"

    if PATHPATTERN in ["IlluminaHiSeq_miRNASeq"]:
        platformTitle ="Illumina HiSeq 2000 RNA Sequencing platform"
    if PATHPATTERN in [ "IlluminaGA_miRNASeq"]:
        platformTitle =" Illumina Genome Analyzer RNA Sequencing platform"

    #change description
    J["description"]=""
    J["dataSubType"]="miRNA isoform expression RNAseq"
    J["label"]= "miRNA isoform expression ("+suffix+")"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") miRNA isoform expression by RNAseq ("+suffix+")"
    J["description"]= J["description"] +"TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") miRNA expression by RNAseq. The miRNA expression profile was measured experimentally using the "+platformTitle+" by the "+ dataProducer +"." + \
        " Level 3 interpreted level data was downloaded from TCGA data coordination center. Download data is in the unit of reads per million mapped reads (RPM). This dataset shows the miRNA transcription estimates in log2 (RPM). For more information see: http://nar.oxfordjournals.org/content/early/2015/08/13/nar.gkv808.full ."

    J["description"] = J["description"] +\
                       "<br><br>In order to more easily view the differential expression between samples, we set the default view to center each miRNA to zero by independently subtracting the mean of each miRNA across the cohort on the fly. Users can view the original non-normalized values by adjusting visualization settings."
    J["description"] = J["description"] +"<br><br>"

    J["wrangling_procedure"]= "Level_3 Data (file names: *.isoform.quantification.txt) download from TCGA DCC, log2(x+1) transformed, and processed at UCSC into Xena repository."

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    J['gdata_tags'] =["transcription","miRNA"]

    #change cgData
    J["name"]="TCGA_"+cancer + "_"+namesuffix
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

    #gene datafile json
    genefile = outDir+cancer+"/"+cgFileName+"_gene"
    oHandle = open(genefile +".json","w")
    J.pop("name")
    J.pop(":probeMap")
    J.pop("description")
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") miRNA gene expression by RNAseq ("+suffix+")"
    J["idSubType"] = "gene"
    J["dataSubType"]="miRNA gene expression RNAseq"
    J["label"]= "miRNA gene expression ("+suffix+")"
    J["wrangling_procedure"]= "Level_3 Data (file names: *.isoform.quantification.txt) download from TCGA DCC, for each sample, all isoform expression for the same miRNA gene are added together, log2(total_RPM +1) transformed, and deposited at UCSC into Xena repository."
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()

    return

def outputProbeMap (outfile, mapping):
    #probeMap
    fout = open(outfile,"w")
    fout.write("#id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
    for gene in mapping:
        for probe in mapping[gene]:
            #probe =hg19:9:139565069-139565091:+
            data = string.split(probe, ":")
            chr = "chr"+data[1]
            start, end = string.split(data[2],'-')
            strand = data[3]
            fout.write(probe+'\t'+gene+'\t'+chr+'\t'+start+'\t'+end+'\t'+strand+'\n')
    fout.close()


def process(dataMatrix, mapping, samples, sample,genes, cancer,infile,flog, valuePOS, LOG2, maxLength):
    # one sample a file
    fin=open(infile,'U')
    fin.readline()
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        hugo = data[1] #miRNA isoform id
        value= data[valuePOS]

        MIMAT = string.split(data[5],',')
        if len(MIMAT)!=2 or MIMAT[0] != "mature":
            continue
        MIMAT = MIMAT[1]
        if MIMAT not in mapping:
            mapping[MIMAT]=[]
        if hugo not in mapping[MIMAT]:
            mapping[MIMAT].append(hugo)

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
                    value = math.log(float(value+1),2)

            x=genes[hugo]
            y=samples[sample]
            dataMatrix[x][y]=value

    fin.close()
    return

def process_percentileRANK (dataMatrix,samples, sample,genes, cancer,infile,flog, valuePOS, maxLength):
    # one sample a file
    fin=open(infile,'U')
    fin.readline()
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        hugo = data[0]
        value= data[valuePOS]
        hugo = string.split(hugo,"|")[0]

        if hugo=="?":
            hugo=data[0]

        if hugo not in genes:
            p=len(genes)
            genes[hugo]=p
            l=[]
            for j in range (0,maxLength):
                l.append("")
            dataMatrix.append(l)
    fin.close()

    #build data from file
    dataDic={}
    fin=open(infile,'U')
    fin.readline()
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        hugo = data[0]
        value= data[valuePOS]
        hugo = string.split(hugo,"|")[0]

        if hugo=="?":
            continue
        try:
            dataDic[hugo]=float(value)
        except:
            dataDic[hugo]="NA"
    fin.close()

    # percentileRANK data
    dataRankDic = percentileRANK (dataDic)

    for hugo in dataRankDic:
        x=genes[hugo]
        y=samples[sample]
        value= dataRankDic[hugo]
        dataMatrix[x][y]=value

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
