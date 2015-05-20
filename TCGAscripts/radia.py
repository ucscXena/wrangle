import sys,string,os
import json,datetime
import math
import inspect
import copy

LEVEL="Level_2"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmpTry/"

#def ucsc_illuminaga_dnaseq_cont_automated_vcf (inDir, outDir, cancer,flog,REALRUN):
#    print cancer, sys._getframe().f_code.co_name
#    distribution = False
#    PATHPATTERN= "IlluminaGA_DNASeq_Cont_automated."
#    PLATFORM = "IlluminaGA"
#    suffix     = "ucsc"
#    namesuffix = "mutation_ucsc_vcf"
#    dataProducer = "University of Californis Santa Cruz GDAC"
#    clean=0
#    radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix, dataProducer,REALRUN,clean, PLATFORM,distribution)

def radiaToXena (inDir, outDir, cancer,flog, PATHPATTERN, suffix, namesuffix,
    dataProducer,REALRUN,clean, PLATFORM,distribution, type):
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
        # go through rootDir file and delete any file WGS in its comment line of DNA_TUMOR
        for dataDir in os.listdir(rootDir):
            badfiles=[]
            for file in os.listdir(rootDir+dataDir):
                fin = open(rootDir+dataDir+"/"+file,'r')
                BAD=0
                while 1:
                    line =fin.readline()
                    if line and line[0]!="#":
                        break
                    if string.find(line,"=WGS") !=-1:
                        BAD=1
                if BAD:
                    badfiles.append(rootDir+dataDir+"/"+file)
            for file in badfiles:
                os.system("rm -f "+ file)

        good=0
        xena = outDir+cancer+"/"+cgFileName
        fout = open (xena,'w')
        fout.write("#"+string.join(["sample","chr","start","end","reference","alt","gene","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")

        for dataDir in os.listdir(rootDir):
            tmpout= "xena_out"
            open(tmpout,'w').close()
            if type =="somatic": #somatic SNP
                os.system("python /inside/home/jzhu/scripts/vcfXenaData/browserDataMelisssa/somaticMutationsForCavm/scripts/runSnpEffOnCohortParseOutput.py -passingSomatic=2 "+ rootDir+dataDir  + " " +tmpout)
            elif type == "germline":
                os.system("python /inside/home/jzhu/scripts/vcfXenaData/browserDataMelisssa/somaticMutationsForCavm/scripts/runSnpEffOnCohortParseOutput.py -passingSomatic=3 "+ rootDir+dataDir  + " " +tmpout)

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

        #get rid of cohort with too many bad samples really stupid UCSC pipeline thing
        r=os.popen("r=$(cut -f 1 "+ xena+ " | sort |uniq| wc -l); echo $r").read()
        if float(string.strip(r))<=10:
            good=0
            print "too many bad samples"

        if good:
            #nonSilentMatrix
            matrixfileout = xena+ "_gene"
            os.system("python xenaToMatrix.py "+ matrixfileout + " "+ xena)
            #nonSilentMatrix json
            nonSilentMatrixJson (matrixfileout+".json", inDir, suffix, cancer, namesuffix+"_gene", dataProducer,PLATFORM, PATHPATTERN)
        else:
            os.system("rm -f "+ xena)
            os.system("rm -f "+ xena+".json")
            os.system("rm -f "+ xena+"_gene_vcf")
            os.system("rm -f "+ xena+"_gene_vcf.json")


    if not os.path.exists(outDir+cancer+"/"+cgFileName):
        return

    #nonSilentMatrix json
    if os.path.exists(outDir+cancer+"/"+cgFileName+"_gene"):
        matrixfileout = outDir+cancer+"/"+cgFileName+ "_gene"
        nonSilentMatrixJson (matrixfileout+".json", inDir, suffix, cancer, namesuffix+"_gene", dataProducer,PLATFORM, PATHPATTERN)

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    J={}
    #stable
    J["cgDataVersion"]=1
    J["dataSubType"]="somatic mutation (SNPs and small INDELs)"
    J["redistribution"]= distribution
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
        J["method"]= "RADIA"
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .vcf file from TCGA DCC, select somatic mutations overlapping with exon region, removed all calls from WGS samples, processed into UCSC Xena muation format, stored in the UCSC Xena repository"

    #change description
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["label"]= "somatic mutation SNPs and small INDELs ("+suffix+" curated vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+" curated vcf)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["label"]= "somatic mutation SNPs and small INDELs ("+suffix+" automated vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+" automated vcf)"
    else:
        J["label"]= "somatic mutation SNPs and small INDELs ("+suffix+" vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") nonsilent somatic mutation ("+suffix+" vcf)"

    J["assembly"]="hg19"
    J["wholeGenome"]= True
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]
    J['owner']="TCGA"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data. Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using the "+ J["method"] +" method."

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
    J["dataSubType"]="somatic non-silent mutation (gene-level)"
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
        J["method"]= "RADIA"
    else:
        J["method"]= ""

    #multiple dirs
    J["url"]=TCGAUtil.remoteBase \
              +string.replace(inDir,TCGAUtil.localBase,"")
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"] ="Download .vcf files from TCGA DCC, select somatic mutations overlapping with exon region, removed all calls from WGS samples, processed into gene by sample matrix of non-silent mutations at UCSC, stored in the UCSC Xena repository"

    #change
    if string.find(PATHPATTERN, "curated")!=-1 :
        J["label"]= "somatic gene-level non-silent mutation ("+suffix+" curated vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation from ("+suffix+" curated vcf)"
    elif string.find(PATHPATTERN, "automated")!=-1 :
        J["label"]= "somatic gene-level non-silent mutation ("+suffix+" automated vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" automated vcf)"
    else:
        J["label"]= "somatic gene-level non-silent mutation ("+suffix+" vcf)"
        J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") gene-level nonsilent somatic mutation ("+suffix+" vcf)"

    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['gdata_tags'] = [dataProducer]
    J['owner']="TCGA"

    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data.  Sequencing data are generated on a "+PLATFORM +" system. The calls are generated at "+dataProducer+" using the "+ J["method"] +" method. <br><br> Red (=1) indicates that a non-silent somatic mutation (nonsense, missense, frame-shif indels, splice site mutations, stop codon readthroughs, change of start codon, inframe indels) was identified in the protein coding region of a gene, or any mutation identified in a non-coding gene. White (=0) indicates that none of the above mutation calls were made in this gene for the specific sample.<br><br>"

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

