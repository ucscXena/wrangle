import sys,string,os
import json,datetime
import math
import inspect
import copy

LEVEL="Level_2"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmptmp/"

#/data/TCGA/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/ucs/gsc/broad.mit.edu/illuminaga_dnaseq_curated/mutations/

def ucsc_illuminaga_dnaseq_cont_vcf (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_Cont."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)
    
def ucsc_illuminaga_dnaseq_cont_automated_vcf (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "IlluminaGA_DNASeq_Cont_automated."
    PLATFORM = "IlluminaGA"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def ucsc_solid_dnaseq_cont_vcf  (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "SOLiD_DNASeq_Cont."
    PLATFORM = "SOLiD"
    suffix     = "ucsc"
    namesuffix = "mutation_ucsc_solid"
    dataProducer = "University of Californis Santa Cruz GDAC"
    clean=1
    radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM)

def radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM):
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
        good=0
        xena = outDir+cancer+"/"+cgFileName 
        fout = open (xena,'w')
        fout.write("#"+string.join(["sample","chr","start","end","gene","reference","alt","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")

        for dataDir in os.listdir(rootDir):
            tmpout= "xena_out"
            open(tmpout,'w').close()
            os.system("python /inside/home/jzhu/scripts/vcfXenaData/browserDataMelisssa/somaticMutationsForCavm/scripts/runSnpEffOnCohortParseOutput.py -passingSomatic=1 "+ rootDir+dataDir  + " " +tmpout)
            
            #remove position not in a gene, bad stupid calls like alt=NA
            fin =open(tmpout,'r')
            for line in fin.readlines():
                if line =="":
                    continue
                elif line[0]=="#":
                    continue
                data = string.split(line,'\t')
                gene = data[4]
                if gene=="":
                    continue
                else:
                    fout.write(line)
                    good=1
            os.system("rm -f "+tmpout)
        fout.close()
        
        if good:
            #nonSilentMatrix
            matrixfileout = xena+ "_gene_vcf"
            os.system("python xenaToMatrix.py "+ matrixfileout + " "+ xena)
            #nonSilentMatrix json
            nonSilentMatrixJson (matrixfileout+".json", inDir, suffix, cancer, namesuffix+"_gene_vcf", dataProducer,PLATFORM, PATHPATTERN)
        else:
            os.system("rm -f "+ xena)

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
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .vcf file from TCGA DCC, processed into UCSC Xena muation format into Xena repository"

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


def nonSilentMatrixJson (jsonFile, inDir, suffix, cancer, namesuffix, dataProducer,PLATFORM, PATHPATTERN):

    oHandle = open(jsonFile,"w")    
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
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .vcf files from TCGA DCC, processed into gene by sample matrix of non-silent mutations at UCSC into Xena repository"
    J["gain"]=10
    J["min"]=-0.1
    J["max"]=0.1

    #change
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+" curated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation from ("+suffix+" curated vcf)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+" automated)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" automated vcf)"
    else:
        J["shortTitle"]= cancer +" gene-level mutation ("+suffix+")"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" vcf)"

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

