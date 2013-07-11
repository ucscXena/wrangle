import string, os, sys, glob
import json,datetime

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

#in the SDRF file, Extract  column will use UUIDs in addition to TCGA barcodes. 
#This can result in an SDRF file that has mixed UUIDs and TCGA barcodes in this column.Namethe 

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/broad.mit.edu/ht_hg-u133a/transcriptome/

PATHPATTERN = "HT_HG-U133A"

def Affy (inDir, outDir, cancer,flog,REALRUN):
    garbage=["tmptmp/"]

    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")
    if os.path.exists( "tmptmp/mageTab/" ):
        os.system("rm -rf tmptmp/mageTab/*")
    else:
        os.system("mkdir tmptmp/mageTab/")

    #multiple files in dir mode
    lastRelease={}
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and  string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
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

        #is tar.gz?, uncompress multiple file mode
        if string.find(file,".tar.gz")!=-1: # and REALRUN:
            os.system("tar -xzf "+inDir+file +" -C tmptmp/") 
            rootDir ="tmptmp/"
            if not REALRUN:
                break

    #print status
    print cancer, __name__

    #mage-tab processing for id mapping
    lastMageTab=0
    for file in os.listdir(inDir):
        if string.find(file,"mage-tab.1")!=-1 and  string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue
        info = string.split(file,".")
        release = int(info [-4])
        if release > lastMageTab:
            lastMageTab =release

    for file in os.listdir(inDir):
        if string.find(file,"mage-tab.1")!=-1 and  string.find(file,".tar.gz")!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue
        info = string.split(file,".")
        release = int(info [-4])
        if release != lastMageTab:
            continue
        if REALRUN:
            os.system("tar -xzf "+inDir+file +" -C tmptmp/mageTab/") 

    mapping={}
    for dirpath, dirnames, filenames in os.walk(rootDir):
        for file in filenames:
            if string.find(file,"sdrf.txt")!=-1:
                mapping= processMageTab(mapping, dirpath+"/"+file)

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

    cgFileName= PATHPATTERN

    #data processing multiple dirs mode
    if REALRUN:
        dataMatrix={}
        samples=[]
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern ="data"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    process(dataMatrix,samples,cancer,infile,mapping, flog)


        outfile = outDir+cancer+"/"+cgFileName
        outputMatrix(dataMatrix, samples, outfile, flog)
    
    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    suffix="AffyU133a"
    J["cgDataVersion"]=1
    J["shortTitle"]="Gene Expression ("+suffix+")"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression ("+suffix+")"
    J[":dataSubType"]="geneExp"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    #J["priority"]= TCGAUtil.browserPriority[cancer]
    J["dataProducer"]= "Broad Institute of MIT and Harvard University cancer genomic characterization center"
    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]=  datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]="tumor"
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["label"]= cancer +" "+J["shortTitle"]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['owner']="TCGA"
    
    #change description
    J["gain"]=1.0
    J["colNormalization"]=True
    J["PLATFORM"]= suffix
    J["wrangling_procedure"]= "Level_3 Data download from TCGA DCC, processed at UCSC into cgData repository"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"\
                      " gene expression data by AffyU133a array.<br><br>"+ \
                      " The gene expression profile was measured experimentally using the Affymetrix HT Human Genome U133a microarray platform by the " + J["dataProducer"] +"."+\
                      " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the gene-level transcription estimates. Data is in log space."+\
                      " Genes are mapped onto the human genome coordinates using UCSC cgData HUGO probeMap."+\
                       "<br><br>In order to more easily view the differential gene expression between samples, we set the default view to center each gene to zero by independently subtracting the mean of the genomic location on the fly. Users can view the original non-normalized values by uncheck the \"Normalize\" option. For more information on how to use the cancer browser, please refer to the help page."
    
    J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
    
    #change cgData
    J["name"]="TCGA_"+cancer+"_exp_u133a"
    name = trackName_fix(J['name'])
    if name ==False:
        message = "bad object name, need fix otherwise break loader, too long "+J["name"]
        print message
        flog.write(message+"\n")
        return
    else:
        J["name"]=name        
    J[":probeMap"]= "hugo"
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
            
    cleanGarbage(garbage)
    return

def processMageTab (dic,file):
    P=7
    
    fin= open(file,'r')
    fin.readline()
    for line in fin.readlines():
        TCGAid = string.split(line,'\t')[0]
        Boardid = string.split(line,'\t')[P]
        if dic.has_key(Boardid) and dic[Boardid]!=TCGAid:
            print  "error",Boardid, dic[Boardid], TCGAid
        dic[Boardid]=TCGAid

    return dic

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf "+dir+"*")
    return

def process(dataMatrix,samples,cancer,infile,mapping, flog):
    # one sample a file
    fin=open(infile,'U')
    line = fin.readline()

    sample = string.split(string.strip(line),"\t")[1]
    if mapping.has_key(sample):
        sample=mapping[sample]
    else:
        fin.close()
        message =  "ERROR sample not in sdrf = "+ sample+ " " +cancer+" "+ __name__
        flog.write(message+"\n")
        print message
        return
    
    if sample in samples:
        fin.close()
        message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__
        flog.write(message+"\n")
        print message
        return

    # Test for barcode or UUID     #throw out all normals and control Analyte
    if sample[0:4]!="TCGA":
        if TCGAUtil.UUID_CELLLINE.has_key(sample):
            print "control cell line ignore", sample
            fin.close()
            return
    else:
        sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
        if sampleTypeCode == False: # likely a uuid
            fin.close()
            return
        elif sampleTypeCode in ["20"]:
            fin.close()
            print "control cell line ignore", sample
            return

    samples.append(sample)
    fin.readline()
    
    for line in fin.readlines():
        hugo,value =string.split(line[:-1],"\t")
        if not dataMatrix.has_key(hugo):
            dataMatrix[hugo]={}
        if value not in ["","null","NULL","Null","NA"]:
            dataMatrix[hugo][sample]=value
        else:
            dataMatrix[hugo][sample]="NA"

    fin.close()
    return 

def  outputMatrix(dataMatrix, samples, outfile, flog):
    fout = open(outfile,"w")
    fout.write("sample\t")
    fout.write(string.join(samples,"\t")+"\n")

    genes = dataMatrix.keys()
    for gene in genes:
        fout.write(gene)
        for sample in samples:
            fout.write("\t"+dataMatrix[gene][sample])
        fout.write("\n")
    fout.close()
