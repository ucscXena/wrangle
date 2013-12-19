import sys,string,os
import json,datetime
import math
import inspect
import copy

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

tmpDir="tmpTry/"

#/inside/home/cline/projects/PanCancer/mutationMatrices/*_cleaned_filtered.maf

def mutationMatrix (inDir, outDir, cancer,flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    PATHPATTERN= "_cleaned_filtered.maf"
    namesuffix = "mutation"
    dataProducer = "TCGA PANCAN AWG"

    garbage=[tmpDir]

    if os.path.exists( tmpDir ):
        os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    #inDir is a single file
    rootDir =""
    remoteDataDirExample =""
    #find the file
    file = inDir

    if string.find(file,PATHPATTERN)!=-1 and string.find(string.upper(file),cancer)!=-1: 
        pass
    else:
        return

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    cgFileName= namesuffix 

    if REALRUN:
        print file, cancer, outDir+cancer+"/"+cgFileName
        fgene=open("/data/TCGA/tcgaDataOneOff/Genomic/PANCAN/genes",'r')
        allGenes=[]
        for gene in string.split(fgene.read(),"\n"):
            allGenes.append(string.strip(gene))
        fgene.close()
        #cut -f 1,9,16 brca_cleaned_filtered.maf |grep $GENE
        process (file, cancer, outDir+cancer+"/"+cgFileName, allGenes)

    oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
    
    J={}
    #stable
    J["cgDataVersion"]=1
    J[":dataSubType"]="somaticMutation"
    J["redistribution"]= True
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J["dataProducer"]= dataProducer
    J["type"]= "genomicMatrix" 
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #multiple dirs
    J["url"]="https://www.synapse.org/#!Synapse:syn1729383"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J[":probeMap"]= "hugo"
    J["shortTitle"]= cancer +" mutation"
    J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation"
    J["gain"]=10
    
    J["label"] = J["shortTitle"] 
    J["sample_type"]=["tumor"]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"

    if cancer!="PANCAN":
        J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
        J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    else:
        J["primary_disease"]="cancer"
        J["anatomical_origin"]=""
        origin =[]
        for value in TCGAUtil.anatomical_origin.values():
            if value =="":
                continue
            if value not in origin:
                origin.extend(value)
        J["tags"]= J["tags"]+origin

    
    J["description"]= "TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") somatic mutation data. Red (=1) indicates that a non-silent somatic mutation (nonsense, missense, frame-shif indels, splice site mutations, stop codon readthroughs) was identified in the protein coding region of a gene, or any mutation identified in a non-coding gene. White (=0) indicates that none of the above mutation calls were made in this gene for the specific sample.<br><br>Somatic mutations calls (even on the same tumor DNA extract) are affected by many factors including library prep, sequencing process, read mapping method, reference genome used, genome annotation, calling algorithms, and ad-hoc pre/postprocessing such as black list genes, target selection regions, and black list samples. This dataset is the best effort made by the TCGA PANCANCER Analysis Working Group. Raw PANCAN mutation data can be downloaded at the url site shown below."
    J["description"] = J["description"] +"<br><br>"
    
    J["wrangling_procedure"]= "TCGA PANCAN strictly filtered maf files (file names: *_cleaned_filtered.maf) download from Synapse, processed into gene by sample matrix at UCSC into cgData repository"

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

def process (file, cancer, outfile, allGenes):
    fin =open(file, 'r')
    fin.readline()
    
    typeDic={
        "Nonsense_Mutation":1,
        "Frame_Shift_Del":2,
        "Frame_Shift_Ins":3,
        "Splice_Site":4,
        "Missense_Mutation":5,
        "Nonstop_Mutation":6,
        "In_Frame_Del":7, ##
        "In_Frame_Ins":8, ##
        "RNA":9,
        "Silent":10 ##
        }
        
    samples=[]
    genes=[]
    dic={}
    for gene in allGenes:
        dic[gene]={}
        
    c=0
    for line in fin.readlines():
        #cut -f 1,9,16 brca_cleaned_filtered.maf |grep $GENE
        c = c+1
        if c % 500 ==0:
            print c
        data =string.split(line[:-1],'\t')
        gene = data[0]
        mtype = typeDic[data[8]]
        sample = data[15]
        if sample not in samples:
            samples.append(sample)
        if gene not in genes:
            genes.append(gene)
        if gene not in dic:
            dic[gene]={}
        if sample not in dic[gene]:
            dic[gene][sample]=mtype
        else:
            if mtype < dic[gene][sample]:
                dic[gene][sample]=mtype
    fin.close()

    print len(genes), len(dic)
    fout=open(outfile,'w')
    fout.write("sample\t"+string.join(samples,"\t")+"\n")
    for gene in dic:
        fout.write(gene)
        for sample in samples:
            try:
                mtype = dic[gene][sample]
                if mtype not in [7,8,10]:
                    fout.write("\t1")
                else:
                    fout.write("\t0")
            except:
                fout.write("\t0")
        fout.write("\n")
    fout.close()
