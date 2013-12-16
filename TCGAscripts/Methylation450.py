import string, os, sys, glob
import json,datetime
import inspect
import copy

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

        
def humanmethylation (inDir, outDir, cancer,flog,PATHPATTERN,BETAOFFSET,REALRUN):
    garbage=["tmptmp/"]
    os.system("rm -rf tmp_*")
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
        allSamples={}
        c=0
        dataMatrix=[]
        tmpSamples={}
        probes={}
        files=[]
        
        oldprobes={}

        GOOD=1

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =PATHPATTERN
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    new_c = process(c, dataMatrix,allSamples,tmpSamples, probes, cancer,infile,flog, BETA_POS,BETAOFFSET,100)
                    #skip
                    if new_c == c:
                        continue
                    c = new_c
                    if (c % 100)==0:
                        tmpout="tmp_"+ str(int(c/100.0))
                        files.append(tmpout)
                        r =outputMatrix(dataMatrix, tmpSamples, probes, oldprobes, tmpout, flog)
                        if r:
                            GOOD=0

                        dataMatrix=[]
                        tmpSamples={}
                        oldprobes=copy.deepcopy(probes)
                        probes ={}

        if (c % 100)!=0:
            tmpout= "tmp_"+ str(int(c/100.0)+1)
            files.append(tmpout)
            r= outputMatrix(dataMatrix, tmpSamples, probes, oldprobes,tmpout, flog)
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
    suffix=PATHPATTERN
    
    J["cgDataVersion"]=1
        
    if suffix =="HumanMethylation450":
        J["shortTitle"]= cancer +" DNA methylation (Methylation450k)"
        namesuffix="hMethyl450" 

    J["label"] = J["shortTitle"] 
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    
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
    
    #change description
    J["gain"]=2.0
    if suffix =="HumanMethylation27":
        J["PLATFORM"]= "Illumina Infinium HumanMethylation27"
    elif suffix =="HumanMethylation450":
        J["PLATFORM"]= "Illumina Infinium HumanMethylation450"
        
    J["wrangling_procedure"]= "Level_3 data download from TCGA DCC, processed at UCSC into cgData repository"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"\
                      " DNA methylation data.<br><br>"+ \
                      " DNA methylation profile was measured experimentally using the "+J["PLATFORM"]+" platform. Beta values were derived at the Johns Hopkins University and University of Southern California TCGA genome characterization center. DNA methylation values, described as beta values, are recorded for each array probe in each sample via BeadStudio software. DNA methylation beta values are continuous variables between 0 and 1, representing the ratio of the intensity of the methylated bead type to the combined locus intensity. Thus higher beta values represent higher level of DNA methylation, i.e. hypermethylation and lower beta values represent lower level of DNA methylation, i.e. hypomethylation.\n\n We observed a bimodal distribution of the beta values from both methylation27 and methylation450 platforms, with two peaks around 0.1 and 0.9 and a relatively flat valley around 0.2-0.8.  The bimodal distribution is far more pronounced and balanced in methylation450 than methylation27 platform. In the methylation27 platform, the lower beta peak is much stronger than the higher beta peak, while the two peaks are of similar hight in the methylation450 platform. During data ingestion to UCSC cgData repository, the beta values were offset by " + str(BETAOFFSET)+" to shift the whole dataset to values between -0.5 to +0.5."
    
    if suffix =="HumanMethylation450":
        J["description"]= J["description"] +"."
        J["description"]= J["description"] +" Microarray probes are mapped onto the human genome coordinates using cgData probeMap derived from GEO GPL13534 record."
    J["description"]=J["description"]+" Here is a <a href=\"http://www.illumina.com/documents/products/appnotes/appnote_dna_methylation_analysis_infinium.pdf\" target=\"_blank\"><u>reference</u></a> to Illumina Infinium BeadChip DNA methylation platform beta value."

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
        os.system("rm -rf "+dir+"*")
    return


def process(c, dataMatrix,allSamples, samples, probes, cancer,infile,flog, BETA_POS, offset, maxLength):
    # one sample a file
    fin=open(infile,'U')    
    line = fin.readline()

    sample = string.split(line[:-1],"\t")[BETA_POS]
    if sample in allSamples or sample in samples:
        fin.close()
        message =  "ERROR duplicated sample = "+ sample+ " " +cancer+" "+ __name__
        flog.write(message+"\n")
        print message
        return c

    # Test for barcode or UUID     #throw out all normals and control Analyte
    if sample[0:4]!="TCGA":
        if TCGAUtil.UUID_CELLLINE.has_key(sample):
            print "control cell line ignore", sample
            fin.close()
            return c
    else:
        sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
        if sampleTypeCode == False: # likely a uuid
            fin.close()
            return c
        elif sampleTypeCode in ["20"]:
            print "control cell line ignore", sample
            fin.close()
            return c

    p=len(samples)
    samples[sample]=p
    allSamples[sample]=""
    c= c+1

    fin.readline()
    
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        probe = data[0]
        value= data[BETA_POS]

        if probe not in probes:
            p=len(probes)
            probes[probe]=p
            l=[]
            for j in range (0,maxLength):
                l.append("")    
            dataMatrix.append(l)
        
        if value not in ["","null","NULL","Null","NA"]:
            value = float(value)+ offset
            x=probes[probe]
            y=samples[sample]
            dataMatrix[x][y]=value

    fin.close()
    return c


def outputMatrix(dataMatrix, samples, probes, oldprobes,outfile, flog):
    #compare probes and oldprobes:
    if oldprobes!={}:
        if len(probes)!=len(oldprobes):
            print "ERROR probes total length is different", len(probes), len(oldprobes), samples
            return 1
        for probe in probes:
            if probes[probe]!=oldprobes[probe]:
                print "ERROR probe order is different", probe
                return 1
            
    fout = open(outfile,"w")
    if oldprobes=={}:
        fout.write("sample")
    for sample in samples:
        fout.write("\t"+sample)
    fout.write("\n")

    for probe in probes:
        if oldprobes=={}:
            fout.write(probe)
        for sample in samples:
            value = dataMatrix[probes[probe]][samples[sample]]
            if value !="":
                value = "%.4f" % (float(value))
                fout.write("\t"+value)
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
    return 0
