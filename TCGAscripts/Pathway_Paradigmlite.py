import string, os, sys
import json,datetime

PATHPATTERN= "Pathway_Paradigm_Lite"
LEVEL="Level_4"
OFFSET =-0.5

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

def Pathway_Paradigmlite (inDir, outDir, cancer,flog,REALRUN):
    if cancer in ["PANCANCER","PANCAN8"]:
        return
    garbage=["tmptmp/"]
    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")

    #figure out the FH date
    if inDir[-1]!="/":
        s = inDir+"/"
    else:
        s=inDir
    FHdate = string.split(s,"/")[-2]

    #single file in dir mode, uncompress to dir
    dataDir =""
    lastDate=""
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue

        if not os.path.exists(inDir +file+".md5"):
            print "file has no matching .md5 throw out", file
            continue

        #file date
        lastDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)
        
        #is tar.gz?, uncompress
        if string.find(file,".tar.gz")!=-1:
            if REALRUN:
                os.system("tar -xzf "+inDir+file +" -C tmptmp/") 
            dataDir ="tmptmp/"+string.replace(file,".tar.gz","")+"/"
            if REALRUN:
                os.system("unzip "+ dataDir+"Evidence_*.zip  -d "+ dataDir)
            break

    #make sure there is data
    if dataDir =="" or (REALRUN and not os.path.exists(dataDir)):
        cleanGarbage(garbage)
        return

    #print status
    print cancer, __name__
    
    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    #data processing single dir mode
    for pattern in ["CNV","Expression","pathlette_IPA-paradigmlite"]:
        file=""
        if REALRUN:
            for file in os.listdir(dataDir):
                if string.find(file,pattern)!=-1:
                    break
        errorTag = cancer +"_"+__name__
        if pattern == "CNV":
            cgFileName= "Paradigm_Evidence_CNV"
            if REALRUN and file!="":
                infile= dataDir+file
                tmpfile="~tmp"
                outfile = outDir+cancer+"/"+cgFileName
                transpose(infile,tmpfile, OFFSET, errorTag,flog)
                process(tmpfile,outfile,cancer)
        if pattern == "Expression":
            cgFileName= "Paradigm_Evidence_Expression"
            if REALRUN and file!="":
                infile= dataDir+file
                tmpfile="~tmp"
                outfile = outDir+cancer+"/"+cgFileName
                transpose(infile,tmpfile, OFFSET, errorTag,flog)
                process(tmpfile,outfile,cancer)
        if pattern == "pathlette_IPA-paradigmlite":
            cgFileName= "Paradigmlite"
            if REALRUN and file!="":
                infile= dataDir+file
                outfile = outDir+cancer+"/"+cgFileName
                process(infile,outfile,cancer)

        oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")
        J={}
        #stable
        if pattern=="CNV":
            suffix="PDMLeCNV"
            J["shortTitle"]="Paradigmlite.evidence.CNV"
            J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") paradigmlite CNV data input"
            J[":dataSubType"]="cna"
            J["gain"]=2.0
        if pattern=="Expression":
            suffix="PDMLeExp"
            J["shortTitle"]="Paradigmlite.evidence.expression"
            J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") paradigmlite gene expression data input"
            J[":dataSubType"]="geneExp"
            J["gain"]=2.0
        if pattern=="pathlette_IPA-paradigmlite":
            suffix="PDML"
            J["shortTitle"]="Paradigmlite"
            J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") paradigmlite inference activity"
            J[":dataSubType"]="PARADIGM.pathlette"                    
            J["gain"]=1.0
                    
        J["cgDataVersion"]=1
        J["redistribution"]= True
        J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
        #J["priority"]= TCGAUtil.browserPriority[cancer]
        J["dataProducer"]= "TCGA FIREHOSE pipeline"
        J["url"]= "http://gdac.broadinstitute.org/runs/analyses__"+FHdate[0:4]+"_"+FHdate[4:6]+"_"+FHdate[6:8]+"/data/"+cancer+"/"+FHdate[0:8]+"/gdac.broadinstitute.org_"+cancer+".Pathway_Paradigm_Lite.Level_4."+FHdate+".0.0.tar.gz"
#        J["version"]= datetime.date.today().isoformat()
        J["version"]=        "2012-06-19"
        J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
                
        #change description
        J["wrangling_procedure"]= "FIREHOSE data download from TCGA DCC, processed at UCSC into cgData repository"
                
        if pattern=="CNV":
            J["description"]= "Broad FireHose automated Paradigmlite run data of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+"). "+\
                              "The dataset shows the input copy number (CNV) data used to derive the PARADIGMLITE inferred gene activities."+\
                              " CNV data were first rank-transformed and fit into a z distribution with final data ranges between 0 to 1 before passing to PARADIGMLITE for analysis. For visualization purpose, these transformed CNV input values were offset by -0.5 to provide a dataset with approximately zero mean and roughly half the values above zero and half the values below zero. Cancer browser uses red/blue coloring for genomic data values above/below zero, providing a visual cue for copy number amplification/deletion."+\
                              " PARADIGM is pathway analysis method to infer patient or sample-specific genetic activities by incorporating curated pathway interactions as well as integrating diverse types of genomic data, such as gene expression and copy number data. PARADIGMLITE, a less computationally intensive version of PARADIGM, integrates different type data on individual genes without incorporating the curated pathway interactions."+\
                              "This dataset shows the input copy number data used for PARADIGMLITE analysis; the PARADIGMLITE inferred gene activities are shown in the accompanied PARADIGMLITE dataset within the same dataset group."+\
                              " Genes are mapped onto the human genome coordinates using cgData HUGO probeMap."+\
                              " Reference to PARADIGM pathway analysis method: PMID 20529912."

        if pattern=="Expression":
            J["description"]= "Broad FireHose automated Paradigmlite run data of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+"). "+\
                              "The dataset shows the input gene expression data used to derive the PARADIGMLITE inferred gene activities. in TCGA "+\
                              " Gene expression data in the form of gene RPKM values, or microarray derived values, were first normalized by the signal from solid normal samples or median/mean normalized per gene when normal samples are not available. The normalized values were then rank-transformed and fit into a z distribution with final data ranges between 0 to 1 before passing to PARADIGMLITE for analysis. For visualization purpose, the transformed gene expression input values were offset by -0.5 to produce a dataset with approximately zero mean and roughly half the values above zero and half the values below zero. Cancer browser uses red/blue coloring for genomic data values above/below zero, providing a visual cue for gene expression over-/under-expression."+ \
                              " PARADIGM is pathway analysis method to infer patient or sample-specific genetic activities by incorporating curated pathway interactions as well as integrating diverse types of genomic data, such as gene expression and copy number data. PARADIGMLITE, a less computationally intensive version of PARADIGM, integrates different type data on individual genes without incorporating the curated pathway interactions."+\
                              " This dataset shows the input gene expression data used for PARADIGMLITE analysis; the PARADIGMLITE inferred gene activities are shown in the accompanied PARADIGMLITE dataset within the same dataset group."+\
                              " Genes are mapped onto the human genome coordinates using cgData HUGO probeMap."+\
                              " Reference to PARADIGM pathway analysis method: PMID 20529912."

        if pattern=="pathlette_IPA-paradigmlite":
            J["description"]= "Broad FireHose automated Paradigmlite run data of TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+"). "+\
                              " Integrated gene activity level inferred using the PARADIGMLITE method."+ \
                              " PARADIGM is pathway analysis method to infer patient or sample-specific genetic activities by incorporating curated pathway interactions as well as integrating diverse types of genomic data, such as gene expression and copy number data. PARADIGMLITE, a less computationally intensive version of PARADIGM, integrates different type data on individual genes without incorporating the curated pathway interactions."+\
                              " This dataset shows the inferred gene activity derived using PARADIGMLITE by integrating gene expression and copy number data; the gene expression and copy number input data are shown in the accompanied evidence datasets within the same dataset group."+\
                              " Genes are mapped onto the human genome coordinates using cgData HUGO probeMap."+\
                              " Reference to PARADIGM pathway analysis method: PMID 20529912."
        J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
        
        #change cgData
        J["name"]="TCGA_"+cancer+"_"+suffix
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

def transpose(infile, outfile, offset, errorTag,flog):
    fin= open(infile,'U')
    lines = fin.readlines()

    ROW=len(lines)
    COL=len(string.split(lines[0],"\t"))

    fin.close()

    matrix =[]
    for i in range (0, COL):
        matrix.append([])
    for i in range (0,COL):
        for j in range (0,ROW):
            matrix[i].append("NA")

    for i in range (0, ROW):
        data = string.split(lines[i][:-1],"\t")
        if len (data)!= COL:
            message= "ERROR: differnt col number in differnet line "+ errorTAG
            print message
            flog.write(message+"\n")
            return 1
        for j in range (0, COL):
            matrix [j][i] =data[j]
    fout= open(outfile,"w")

    fout.write(string.join(matrix[0],"\t")+"\n")
    for i in range (1,COL):
        fout.write(matrix[i][0])
        for j in range (1,ROW):
            value = matrix[i][j]
            if value != "NA":
                value= float(value)+offset
                fout.write("\t"+str(value))
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
    return 0

def process(infile, outfile,cancer):
    if cancer in ["LAML"]:
        cancerCode="03"
    else:
        cancerCode="01"
        
    fin=open(infile,'r')
    fout=open(outfile,'w')

    line =fin.readline()[:-1]
    ids= string.split(line,'\t')
    fout.write(ids[0])
    for id in ids[1:]:
        fout.write("\t"+id+ "-"+cancerCode)
    fout.write("\n")

    for line in fin.readlines():
        fout.write(line)
    fin.close()
    fout.close()
    
