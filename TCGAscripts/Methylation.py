import string, os, sys, glob
import json,datetime
import inspect

LEVEL="Level_3"
BETAOFFSET =-0.5
BETA_POS=1

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/jhu-usc.edu/humanmethylation450/methylation/
#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/jhu-usc.edu/humanmethylation27/methylation/

def humanmethylation450 (inDir, outDir, cancer,flog,REALRUN):
    print cancer, inspect.stack()[0][3]     
    PATHPATTERN= "HumanMethylation450"
    humanmethylation (inDir, outDir, cancer,flog, PATHPATTERN,BETAOFFSET,REALRUN)
    return

def humanmethylation27 (inDir, outDir, cancer,flog,REALRUN):
    print cancer, inspect.stack()[0][3] 
    PATHPATTERN= "HumanMethylation27"
    humanmethylation (inDir, outDir, cancer,flog, PATHPATTERN,BETAOFFSET,REALRUN)
    return
        
def humanmethylation (inDir, outDir, cancer,flog,PATHPATTERN,BETAOFFSET,REALRUN):
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
    if lastDate == None:
        cleanGarbage(garbage)
        print "ERROR expect data, but no data in dirpath", rootDir, cancer, __name__
        return
    
    if REALRUN and (rootDir =="" or not os.path.exists(rootDir)):
        cleanGarbage(garbage)
        print "ERROR expect data, but wrong dirpath", rootDir, cancer, __name__
        return
    
    #print status
    print cancer, __name__
    
    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= PATHPATTERN

    #data processing multiple dirs mode
    if REALRUN:
        samples={}
        probes={}
        total =0.0
        count=0
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =PATHPATTERN
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    total, count = betaMean(total, count,samples, probes, cancer,infile,flog)
        average = total/count
        print average,BETAOFFSET
    
        dataMatrix=[]
        for i in range(0,len(probes)):
            l=[]
            for j in range (0,len(samples)):
                l.append("NA")
            dataMatrix.append(l)

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =PATHPATTERN
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    process(dataMatrix,samples,probes, cancer,infile,flog,BETAOFFSET)

        outfile = outDir+cancer+"/"+cgFileName
        outputMatrix(dataMatrix, samples, probes,outfile, flog)
    
    if not REALRUN:
        iHandle = open(outDir+cancer+"/"+cgFileName+".json","r")
        J = json.loads(iHandle.read())
        iHandle.close()
        average = J["mean"]

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    suffix=PATHPATTERN
    J["cgDataVersion"]=1
    if suffix =="HumanMethylation27":
        J["shortTitle"]="DNA Methylation (Methylation27)"
        namesuffix="hMethyl27"
    if suffix =="HumanMethylation450":
        J["shortTitle"]="DNA Methylation (Methylation450)"
        namesuffix="hMethyl450" 
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") DNA methylation ("+suffix+")"
    J[":dataSubType"]="DNAMethylation"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    #J["priority"]= TCGAUtil.browserPriority[cancer]            
    J["dataProducer"]= "Johns Hopkins University and University of Southern California TCGA cancer genomic characterization center"
    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["valOffset"] = BETAOFFSET
    J["mean"]= average
    #change description
    J["gain"]=2.0
    if suffix =="HumanMethylation27":
        J["PLATFORM"]= "Illumina Infinium HumanMethylation27"
    elif suffix =="HumanMethylation450":
        J["PLATFORM"]= "Illumina Infinium HumanMethylation450"
        
    J["wrangling_procedure"]= "Level_3 data download from TCGA DCC, processed at UCSC into cgData repository"

    J["description"]= "The dataset shows TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"\
                      " DNA methylation data."+ \
                      " DNA methylation profile was measured experimentally using the "+J["PLATFORM"]+" platform. Beta values were derived at the Johns Hopkins University and University of Southern California TCGA genome characterization center. DNA methylation values, described as beta values, are recorded for each array probe in each sample via BeadStudio software. DNA methylation beta values are continuous variables between 0 and 1, representing the ratio of the intensity of the methylated bead type to the combined locus intensity. Thus higher beta values represent higher level of DNA methylation, i.e. hypermethylation and lower beta values represent lower level of DNA methylation, i.e. hypomethylation.\n\n We observed a bimodal distribution of the beta values from both methylation27 and methylation450 platforms, with two peaks around 0.1 and 0.9 and a relatively flat valley around 0.2-0.8.  The bimodal distribution is far more pronounced and balanced in methylation450 than methylation27 platform. In the methylation27 platform, the lower beta peak is much stronger than the higher beta peak, while the two peaks are of similar hight in the methylation450 platform. During data ingestion to UCSC cgData repository, the beta values were offset by " + str(BETAOFFSET)+" to shift the whole dataset to values between -0.5 to +0.5. The average of the unshifted beta values of this dataset is "+str(average)
    if suffix =="HumanMethylation27":
        J["description"]= J["description"] +", thus much of the heatmap appears hypomethylated (blue)."
        J["description"]= J["description"] +" Microarray probes are mapped onto the human genome coordinates using cgData probeMap derived from GEO GPL8490 and GPL13534 records."
    if suffix =="HumanMethylation450":
        J["description"]= J["description"] +"."
        J["description"]= J["description"] +" Microarray probes are mapped onto the human genome coordinates using cgData probeMap derived from GEO GPL13534 record."
    J["description"]=J["description"]+" Here is a <a href=\"http://www.illumina.com/documents/products/appnotes/appnote_dna_methylation_analysis_infinium.pdf\" target=\"_blank\"><u>reference</u></a> to Illumina Infinium BeadChip DNA methylation platform beta value."

    J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
    
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
    if suffix =="HumanMethylation27":
        J[":probeMap"]= "illuminaMethyl27_GPL8490_GPL13534"
    if suffix =="HumanMethylation450":
        J[":probeMap"]= "illuminaHumanMethylation450_GPL13534"
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
            
    cleanGarbage(garbage)
    return

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf dir")
    return

def betaMean(total, count,samples, probes,cancer,infile,flog):
    # one sample a file
    fin=open(infile,'r')
    line = fin.readline()
    sample = string.split(string.strip(line),"\t")[BETA_POS]
    if sample in samples:
        fin.close()
        message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__
        flog.write(message+"\n")
        print message
        return total, count


    # Test for barcode or UUID     #throw out all normals and control Analyte
    if sample[0:4]!="TCGA":
        print sample
        if TCGAUtil.UUID_CELLLINE.has_key(sample):
            print "control cell line ignore", sample
            fin.close()
            return total,count
    else:
        sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
        if sampleTypeCode == False: # likely a uuid
            fin.close()
            return
        elif sampleTypeCode in ["20"]:
            print "control cell line ignore", sample
            fin.close()
            return total,count

    p=len(samples)
    samples[sample]=p
    
    fin.readline()

    for line in fin.readlines():
        probe,value =string.split(line[:-1],"\t")[0:BETA_POS+1]
        if probe not in probes:
            p = len(probes)
            probes[probe]= p
        if value in ["","null","NULL","Null","NA"]:
            continue
        value = float(value)
        total = total +value
        count = count +1
    fin.close()
    return total, count

def process(dataMatrix,samples, probes,cancer,infile,flog, offset):
    # one sample a file
    fin=open(infile,'r')
    line = fin.readline()
    sample = string.split(string.strip(line),"\t")[BETA_POS]

    if not samples.has_key(sample):
        return
    
    sampleP= samples[sample]
    
    fin.readline()

    for line in fin.readlines():
        probe,value =string.split(line[:-1],"\t")[0:BETA_POS+1]
        probeP = probes[probe]
        if value not in ["","null","NULL","Null","NA"]:
            value = float(value)
            dataMatrix[probeP][sampleP]=value+offset
    fin.close()
    return 

def  outputMatrix(dataMatrix, samples, probes, outfile, flog):
    fout = open(outfile,"w")

    samplesStr=[]
    for i in range (0, len(samples)):
        samplesStr.append("NA")
    for sample in samples:
        p= samples[sample]
        samplesStr[p]=sample
    fout.write("sample\t")
    fout.write(string.join(samplesStr,"\t")+"\n")

    probesStr=[]
    for i in range (0, len(probes)):
        probesStr.append("NA")
    for probe in probes:
        p= probes[probe]
        probesStr[p]=probe

    for i in range (0, len(probesStr)):
        fout.write(probesStr[i])
        for j in range (0, len(samplesStr)):
            if dataMatrix[i][j]!="NA":
                fout.write("\t"+"%.4f" % (dataMatrix[i][j]))
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
