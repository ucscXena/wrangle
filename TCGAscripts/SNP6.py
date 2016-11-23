import string, os, sys, glob
import json,datetime

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

garbage=["tmptmp/"]

#in the SDRF file, Extract  column will use UUIDs in addition to TCGA barcodes. 
#This can result in an SDRF file that has mixed UUIDs and TCGA barcodes in this column.Namethe 

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/broad.mit.edu/genome_wide_snp_6/snp/
PATHPATTERN = "Genome_Wide_SNP_6"

def SNP6 (inDir, outDir, cancer,flog, REALRUN):
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
                mapping= processMageTab(mapping, dirpath+"/"+file,7)

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

    outfile =  outDir+cancer+"/"+"SNP6_genomicSegment"
    gMoutput = outDir+cancer+"/"+"SNP6"
    if REALRUN:
        fout=open(outfile,'w')
        samples=[]
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern =".hg19.seg."
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    process(samples,cancer,infile,flog,mapping, fout)
        fout.close()
        if samples!=[]:
            # segToMatrix
            #os.system("python seg2matrix/segToMatrix.py "+outfile +" seg2matrix/refGene_hg19 "+ gMoutput)
            os.system("python seg2matrix/segToMatrixGalaxy.py "+outfile +" seg2matrix/refGene_hg19 "+ gMoutput+".matrix "+gMoutput+".probeMap 0")
        else:
            os.system("rm -f "+outfile)
            os.system("rm -f "+outfile+".json")
            os.system("rm -f "+gMoutput+".matrix")
            os.system("rm -f "+gMoutput+".probeMap")
            os.system("rm -f "+gMoutput+".matrix.json")
            os.system("rm -f "+gMoutput+".probeMap.json")
            
    noCNV=0
    if os.path.exists(outfile):
        oHandle = open(outfile+".json","w")
        makeJSON(oHandle,cancer,lastMageTab,inDir,noCNV,"genomicSegment")
    if os.path.exists(gMoutput+".matrix"):
        oHandle = open(gMoutput+".matrix.json","w")
        makeJSON(oHandle,cancer,lastMageTab,inDir,noCNV,"genomicMatrix")
        #probeMap
        oHandle = open(gMoutput+".probeMap.json","w")
        makeProbeJSON(oHandle,cancer,noCNV)
    


    outfile =  outDir+cancer+"/"+"SNP6_nocnv_genomicSegment"
    gMoutput = outDir+cancer+"/"+"SNP6_nocnv"
    if REALRUN:
        fout=open(outfile,'w')
        samples=[]
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                pattern ="nocnv_hg19.seg."
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    process(samples,cancer,infile,flog,mapping, fout)
        fout.close()
        if samples!=[]:
            # segToMatrix
            #os.system("python seg2matrix/segToMatrix.py "+outfile +" seg2matrix/refGene_hg19 "+ gMoutput)
            os.system("python seg2matrix/segToMatrixGalaxy.py "+outfile +" seg2matrix/refGene_hg19 "+ gMoutput+".matrix "+gMoutput+".probeMap 0")
        else:
            os.system("rm -f "+outfile)
            os.system("rm -f "+outfile+".json")
            os.system("rm -f "+gMoutput+".matrix")
            os.system("rm -f "+gMoutput+".probeMap")
            os.system("rm -f "+gMoutput+".matrix.json")
            os.system("rm -f "+gMoutput+".probeMap.json")
    noCNV=1
    if os.path.exists(outfile):
        oHandle = open(outfile+".json","w")
        makeJSON(oHandle,cancer,lastMageTab,inDir,noCNV,"genomicSegment")
    if os.path.exists(gMoutput+".matrix"):
        oHandle = open(gMoutput+".matrix.json","w")
        makeJSON(oHandle,cancer,lastMageTab,inDir,noCNV,"genomicMatrix")
        #probeMap
        oHandle = open(gMoutput+".probeMap.json","w")
        makeProbeJSON(oHandle,cancer,noCNV)
            
    cleanGarbage(garbage)
    return

def makeProbeJSON(oHandle,cancer,noCNV):
    J={}
    #change cgData
    if noCNV:
        mName="TCGA_"+cancer+"_GSNP6noCNV"
    else:
        mName="TCGA_"+cancer+"_GSNP6raw"
            
    mName = trackName_fix(mName)
    if mName ==False:
        message = "bad object name, need fix otherwise break loader, too long "+mName
        print message
        flog.write(message+"\n")
        return

    J["name"]=mName+"_probeMap"
    J["type"]="probeMap"
    J["version"]= datetime.date.today().isoformat()
    J["assembly"]="hg19"
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
    
def makeJSON(oHandle,cancer,lastMageTab,inDir,noCNV,type):    
    J={}
    #stable
    if noCNV:
        J["label"]= "copy number (delete germline cnv)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") segmented copy number (delete germline cnv)"
    else:
        J["label"]= "copy number"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") segmented copy number"
    J["dataSubType"]="copy number"
    J["dataProducer"]= "Broad Institute of MIT and Harvard University cancer genomic characterization center"
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"
    J["unit"]= "log2(tumor/normal)"
    #change description
    J["PLATFORM"]= "SNP6"
    J["wrangling_procedure"]= "Level_3 data download from TCGA DCC, processed at UCSC into cgData repository"

    if noCNV:
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") segmented copy number variation profile after removing common germline copy number variation.<br><br>"
    else:
        J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") segmented copy number variation profile.<br><br>"
    J["description"]=J["description"]+" Copy number profile was measured experimentally using the Affymetrix Genome-Wide Human SNP Array 6.0 platform at the Broad TCGA genome characterization center. Raw copy numbers were estimated at each of the SNP and copy-number markers. <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/15475419\" target=\"_blank\"><u>Circular binary segmentation</u></a> was then used to segment the copy number data. Segments are mapped to hg19 genome assembly at Broad."
    if noCNV:
        J["description"] = J["description"]+" A fixed set of common germline cnv probes were removed prior to segmentation."
    J["description"]=J["description"]+" Reference to the algorithm used by Broad to produce the dataset: <a href =\"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"+string.lower(cancer)+"/cgcc/broad.mit.edu/genome_wide_snp_6/snp/broad.mit.edu_"+cancer+".Genome_Wide_SNP_6.mage-tab.1."+ str(lastMageTab)+".0/DESCRIPTION.txt\" target=\"_blank\"> <u>DCC description</u></a> and <a href=\"http://www.nature.com/nature/journal/vaop/ncurrent/suppinfo/nature07385.html\" target=\"_blank\"><u>nature 2008</u></a> ."

    J["description"] = J["description"] +"<br><br>"
                
    #change cgData
    if type=="genomicSegment":
        if noCNV:
            J["name"]="TCGA_"+cancer+"_GSNP6noCNV_gSeg"
        else:
            J["name"]="TCGA_"+cancer+"_GSNP6raw_gSeg"

    if type=="genomicMatrix":
        if noCNV:
            J["name"]="TCGA_"+cancer+"_GSNP6noCNV"
        else:
            J["name"]="TCGA_"+cancer+"_GSNP6raw"
            
    name = trackName_fix(J['name'])
    if name ==False:
        message = "bad object name, need fix otherwise break loader, too long "+J["name"]
        print message
        flog.write(message+"\n")
        return
    else:
        J["name"]=name        

    J["type"]= type
    if type=="genomicSegment":
        J["assembly"]= "hg19"
    if type=="genomicMatrix":
        J[":probeMap"]=J['name']+"_probeMap"
        if noCNV:
            J[":genomicSegment"] =trackName_fix("TCGA_"+cancer+"_GSNP6noCNV_gSeg")
        else:
            J[":genomicSegment"] =trackName_fix("TCGA_"+cancer+"_GSNP6raw_gSeg")

        
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
                
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()

def processMageTab (dic,file,P):
    fin= open(file,'r')
    fin.readline()
    for line in fin.readlines():
        TCGAid = string.split(line,'\t')[0]
        Boardid = string.split(line,'\t')[P]
        if dic.has_key(Boardid) and dic[Boardid]!=TCGAid:
            print  "error",Boardid, dic[Boardid], TCGAid
            print string.split(line,'\t')[0:P+1]
        dic[Boardid]=TCGAid
    return dic

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf "+ dir+"*")
    return

def process(samples,cancer,infile,flog, mapping,fout):
    # one sample a file
    fin=open(infile,'U')
    line = fin.readline()

    line = fin.readline()
    sample = string.split(string.strip(line),"\t")[0]
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
        if TCGAUtil.UUID_NORMAL_CELLLINE.has_key(sample):
            fin.close()
            return
            
    else:
        sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
        if sampleTypeCode == False: # likely a uuid
            fin.close()
            return
        elif sampleTypeCode in ["10","11","20"]:
            #check TCGAUtil for codes
            fin.close()
            return

        
    samples.append(sample)
    fin.close()

    fin=open(infile,'U')
    line = fin.readline()
    fout.write("sample"+"\t"+"chr"+"\t"+"start"+"\t"+"end"+"\t"+"value"+"\n")
    for line in fin.readlines():
        sample,chr,start, end, numMark, segMean = string.split(line[:-1],"\t")
        sample=mapping[sample]
        if chr=="23":
            chr="chrX"
        elif chr=="24":
            chr="chrY"
        elif chr=="M":
            continue
        else:
            chr="chr"+chr
        start = str(int(float(start)))
        end = str(int(float(end)))
        fout.write(sample+"\t"+chr+"\t"+start+"\t"+end+"\t"+segMean+"\n")
