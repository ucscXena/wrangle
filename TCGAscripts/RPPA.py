import sys,string,os
import json,datetime
import math
import inspect
import copy

LEVEL="Level_3"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmpTry/"

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/cgcc/mdanderson.org/mda_rppa_core/protein_exp/

def RPPA (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name

    PATHPATTERN = "MDA_RPPA_Core"    
    dataProducer= "MD Anderson Cancer Center TCGA proteome characterization center"
    
    garbage=[tmpDir]

    if os.path.exists( tmpDir ):
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
        if string.find(file,".tar.gz")!=-1 and REALRUN :
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

    cgFileName= "RPPA"
    
    #data processing multiple dirs mode
    if REALRUN:
        aliquote_dic =TCGAUtil.uuid_Aliquot_all()
        dataMatrix={}
        allSamples=[]
        probes=[]
        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample =""
                pattern ="protein_expression"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    sample = string.split(file,".")[5]
                if sample =="":
                    continue            
                # Test for barcode or UUID     #throw out all normals and control Analyte
                if sample[0:4]!="TCGA":
                    if aliquote_dic.has_key(string.lower(sample)):
                        if TCGAUtil.UUID_CELLLINE.has_key(sample):
                            print "control cell line ignore", sample
                            continue
                    else:
                        print "unknow id:", sample
                        continue
                else:
                    sampleTypeCode = TCGAUtil.barcode_SampleType(sample)
                    if sampleTypeCode == False: # likely a uuid
                        continue
                    elif sampleTypeCode in ["20"]:
                        print "control cell line ignore", sample
                        continue
                if sample not in allSamples:
                    allSamples.append(sample)

        for dataDir in os.listdir(rootDir):
            for file in os.listdir(rootDir+dataDir):
                sample =""
                pattern ="protein_expression"
                if string.find(file,pattern)!=-1:
                    infile = rootDir+dataDir+"/"+file
                    sample = string.split(file,".")[5]
                if sample =="":
                    continue
                if sample not in allSamples:
                    continue
                valuePOS=1
                process(dataMatrix,allSamples,sample, probes, cancer,infile,flog, valuePOS)

    
        outfile = outDir+cancer+"/"+cgFileName
        outputMatrix(dataMatrix, allSamples, probes, outfile)

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    J["dataSubType"]="protein expression RPPA"
    J["redistribution"]= True
    J["dataProducer"]= dataProducer
    J["colNormalization"]=True
    J["PLATFORM"]= "M.D. Anderson Reverse Phase Protein Array Core platform"
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "Xena TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"]= "Level_3 Data (file names: *.protein_expression.*) download from TCGA DCC, and processed at UCSC into Xena repository"
    #J["label"]= "protein expression RPPA"
    J["label"]= "RPPA"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") reverse phase protein array"
    
    J[":probeMap"]= "md_anderson_antibodies"

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['owner']="TCGA"
    J["tags"]=["cancer"]+ TCGAUtil.tags[cancer]
    J["unit"]="normalized RPPA value"
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") protein expression by reverse phase protein array (RPPA).<br><br> The data was generated and processed at the MD Anderson Cancer Center TCGA proteome characterization center RPPA core. Level 3 interpreted level data was downloaded from TCGA data coordination center.<br><br>"
    
    J["description"] = J["description"] + "Data normalization from the MDACC RPPA core: <a href=\"http://bioinformatics.mdanderson.org/main/TCPA:Overview\" target=\"_blank\"><u>under section How are the RPPA data processed</u></a>.<br>"
    
    #change cgData
    J["name"]="TCGA_"+cancer+"_RPPA"
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

def process(dataMatrix,samples, sample,genes, cancer,infile,flog, valuePOS):
    # one sample a file
    fin=open(infile,'U')    
    fin.readline()
    fin.readline()

    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        antibody = data[0]
        value= data[valuePOS]

        if antibody not in genes:
            genes.append(antibody)
            dataMatrix[antibody]={}
            for row in samples:
                dataMatrix[antibody][row]=""

        if value not in ["","null","NULL","Null","NA"]:
            dataMatrix[antibody][sample]= value
    fin.close()
    return 

def outputMatrix(dataMatrix, samples, genes, outfile):
    fout = open(outfile,"w")
    fout.write("sample")
    for sample in samples:
        fout.write("\t"+sample)
    fout.write("\n")

    for gene in genes:
        fout.write(gene)
        for sample in samples:
            value = dataMatrix[gene][sample]
            if value !="":
                fout.write("\t"+str(value))
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
    return 0
