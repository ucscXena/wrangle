import string, os, sys
import json,datetime
import csv
import filecmp

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from CGDataLib import *
import  TCGAUtil

def genomic (dir,outDir, cancer,flog,REALRUN):
    if cancer not in ["LUNG","COADREAD","GBMLGG"]:
        return

    print cancer, sys._getframe().f_code.co_name

    if cancer =="LUNG":
        c1 ="LUSC"
        c2 ="LUAD"
    if cancer =="COADREAD":
        c1 ="COAD"
        c2 ="READ"
    if cancer =="GBMLGG":
        c1 ="LGG"
        c2 ="GBM"

    for file in [
        "mutation",
        "HiSeqV2",
        "HiSeqV2_exon",
        "AgilentG4502A_07_3",
        "Gistic2_CopyNumber_Gistic2_all_data_by_genes",
        "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
        "HumanMethylation27",
        "HumanMethylation450",
        "RPPA_RBN",
        ""
        ]:

        if file =="mutation" and cancer =="COADREAD": #PANCAN AWG mutation
            continue

        type= "matrix"

        process (outDir, cancer, c1, c2, file, REALRUN,type)
    return


def process (outDir, cancer, c1, c2, file, REALRUN,type):
    if file =="":
        return

    c1dir = outDir+c1+"/"
    c2dir = outDir+c2+"/"

    c1file = c1dir+file
    c2file = c2dir+file

    if ( not os.path.exists(c1file)  or  not os.path.exists(c2file)):
        return

    if REALRUN:
        if type =="mutation":
            print "todo"
        elif type =="SNP6_genomicSegment":
            os.system("cat "+c1file+" "+c2file +" > "+outDir+cancer+"/"+file)
            # segToMatrix
            outfile=outDir+cancer+"/"+file
            gMoutput = outDir+cancer+"/"+string.replace(file,"_genomicSegment","")
            os.system("python seg2matrix/segToMatrix.py "+outfile +" seg2matrix/refGene_hg18 "+ gMoutput)
        else: 
            os.system("python mergeGenomicMatrixFiles.py "+outDir+cancer+"/"+file + " "+c1file +" "+ c2file)

    J={}
    iHandle =open(c1file+".json","r")
    Jinput= json.loads(iHandle.read())
    print c1file
    J["dataSubType"]=  Jinput["dataSubType"]
    J["version"]=  Jinput["version"]
    J["type"]=  Jinput["type"]
    s= Jinput["label"]
    s = string.replace(s,c1,cancer)
    J["label"]= s
    if J["type"] not in ["genomicSegment","mutationVector","clinicalMatrix"]:
        J[":probeMap"]=Jinput[":probeMap"]
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"

    if Jinput.has_key("gdata_tags"):
        J["gdata_tags"]= Jinput["gdata_tags"]
    if Jinput.has_key("min"):
        J["min"]= Jinput["min"]
    if Jinput.has_key("max"):
        J["max"]= Jinput["max"]
    if Jinput.has_key("colNormalization"):
        J["colNormalization"]= Jinput["colNormalization"]
        
    s= Jinput ["longTitle"]
    s = string.replace(s,"("+c1+")","("+cancer+")")
    s = string.replace(s,TCGAUtil.cancerOfficial[c1],TCGAUtil.cancerOfficial[cancer])
    J["longTitle"]=s

    s=Jinput["description"]
    s = string.replace(s,"("+c1+")","("+cancer+")")
    s = string.replace(s,TCGAUtil.cancerOfficial[c1],TCGAUtil.cancerOfficial[cancer])
    J["description"]= "The dataset is combined from TCGA "+ TCGAUtil.cancerOfficial[c1]+" and "+ TCGAUtil.cancerOfficial[c2]+" datasets. "+ s
    J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"

    s = Jinput["name"]
    s = string.replace(s,c1,cancer)
    J["name"]=s
    output = open(outDir+cancer+"/"+file+".json","w")
    output.write(json.dumps(J,indent=-1))
    output.close()

    if type =="SNP6_genomicSegment":
        J={}
        gMoutput = outDir+cancer+"/"+string.replace(file,"_genomicSegment","")
        cfile = gMoutput+".matrix.json"
        iHandle =open(cfile,"r")
        Jinput= json.loads(iHandle.read())
        J["dataSubType"]=  Jinput["dataSubType"]
        J["version"]=  Jinput["version"]
        J["type"]=  Jinput["type"]
        J["cgDataVersion"]=     Jinput["cgDataVersion"]
        s= Jinput["shortTitle"]
        s = string.replace(s,c1,cancer)
        J["shortTitle"]= s
        J["label"] = J["shortTitle"] 
        J[":probeMap"]=string.replace(Jinput[":probeMap"],c1,cancer)
        
        if Jinput.has_key("gain"):
            J["gain"]= Jinput["gain"]
        if Jinput.has_key("min"):
            J["min"]= Jinput["min"]
        if Jinput.has_key("max"):
            J["max"]= Jinput["max"]
        
        s= Jinput ["longTitle"]
        s = string.replace(s,c1,cancer)
        s = string.replace(s,TCGAUtil.cancerOfficial[c1],TCGAUtil.cancerOfficial[cancer])
        J["longTitle"]=s
        J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"

        s = Jinput["name"]
        s = string.replace(s,c1,cancer)
        J["name"]=s
        output = open(outDir+cancer+"/"+gMoutput+".matrix.json","w")
        output.write(json.dumps(J,indent=-1))
        output.close()

        J={}
        cfile = c1dir+gMoutput+".probeMap.json"
        iHandle =open(cfile+".json","r")
        Jinput= json.loads(iHandle.read())
        J=Jinput
        J["name"]=string.replace(J["name"],c1, cancer)
        output = open(outDir+cancer+"/"+gMoutput+".probeMap.json","w")
        output.write(json.dumps(J,indent=-1))
        output.close()
    return
