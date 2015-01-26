import string, os, sys, glob
import json,datetime
import inspect

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/

def AgilentG4502A_07_3 (inDir, outDir, cancer,flog,REALRUN):
    print cancer, inspect.stack()[0][3]
    PATHPATTERN= "AgilentG4502A_07_3"
    AgilentG4502A (inDir, outDir, cancer,flog, PATHPATTERN,REALRUN)
    return

def AgilentG4502A_07_2 (inDir, outDir, cancer,flog,REALRUN):
    print cancer, inspect.stack()[0][3]
    PATHPATTERN= "AgilentG4502A_07_2"
    AgilentG4502A (inDir, outDir, cancer,flog, PATHPATTERN,REALRUN)
    return

def AgilentG4502A_07_1 (inDir, outDir, cancer,flog,REALRUN):
    print cancer, inspect.stack()[0][3]
    PATHPATTERN= "AgilentG4502A_07_1"
    AgilentG4502A (inDir, outDir, cancer,flog, PATHPATTERN,REALRUN)
    return
        
def AgilentG4502A (inDir, outDir, cancer,flog,PATHPATTERN,REALRUN):
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
                    process(dataMatrix,samples,cancer,infile,flog)


        outfile = outDir+cancer+"/"+cgFileName
        outputMatrix(dataMatrix, samples, outfile, flog)
    
    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    suffix=PATHPATTERN
    J["cgDataVersion"]=1
    J["shortTitle"]=cancer +" gene expression ("+suffix+")"
    J["label"] = J["shortTitle"] 
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene expression ("+suffix+" array)"
    J[":dataSubType"]="gene expression array"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    #J["priority"]= TCGAUtil.browserPriority[cancer]
    J["dataProducer"]= "University of North Carolina TCGA genome characterization center"
    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]=  datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    J['gdata_tags'] =["transcription"]
    
    #change description
    J["gain"]=1.0
    J["min"]=-1.0
    J["max"]=1.0
    J["colNormalization"]=True    
    J["PLATFORM"]= suffix
    J["wrangling_procedure"]= "Level_3 Data download from TCGA DCC, processed at UCSC into cgData repository"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"\
                      " gene expression data by agilent array.<br><br>"+ \
                      " The gene expression profile was measured experimentally using Agilent 244K custom gene expression "+string.upper(PATHPATTERN[7:])+" microarrays by the University of North Carolina TCGA genomic characterization center. "+\
                      " Level 3 interpreted level data was downloaded from TCGA data coordination center. This dataset shows the gene level transcription estimates (level 3 data), as in log2 lowess normalized ratio of sample signal to reference signal (cy5/cy3) collapsed by gene."+\
                      " Genes are mapped onto the human genome coordinates using UCSC cgData HUGO probeMap."+\
                      "<br><br>In order to more easily view the differential gene expression between samples, we set the default view to be gene-level normalized by independently subtracting the mean of the genomic location on the fly. Users can view the original non-normalized values by uncheck the \"Normalize\" option. For more information on how to use the cancer browser, please refer to the help page."
    
    J["description"] = J["description"] +"<br><br>"
    
    #change cgData
    J["name"]="TCGA_"+cancer+"_"+string.replace(suffix,"Agilent","")
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

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf "+dir+"*")
    return

def process(dataMatrix,samples, cancer,infile,flog):
    # one sample a file
    fin=open(infile,'U')
    line = fin.readline()
    sample = string.split(string.strip(line),"\t")[1]
    if sample in samples:
        fin.close()
        message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__
        flog.write(message+"\n")
        print message
        return

    # Test for barcode or UUID     #throw out all normals and control Analyte
    if sample[0:4]!="TCGA":
        print sample
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
