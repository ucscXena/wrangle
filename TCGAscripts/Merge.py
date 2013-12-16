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
    if cancer not in ["LUNG","COADREAD"]:
        return

    print cancer, sys._getframe().f_code.co_name

    if cancer =="LUNG":
        c1 ="LUSC"
        c2 ="LUAD"
    if cancer =="COADREAD":
        c1 ="COAD"
        c2 ="READ"

    for file in [
        "HiSeqV2",
        "HiSeqV2_exon",
        "AgilentG4502A_07_3",
        "Gistic2_CopyNumber_Gistic2_all_data_by_genes",
        "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
        "Gistic2_CopyNumber_Gistic2_focal_data_by_genes",
        "HumanMethylation27",
        "HumanMethylation450",
        "RPPA",
        "RPPA_RBN",
        #"SNP6_genomicSegment",
        #"SNP6_nocnv_genomicSegment",
        ""
        ]:
        
        type=""
        if file =="RPPA":
            type="RPPA"
        if file=="SNP6_genomicSegment" or file =="SNP6_nocnv_genomicSegment":
            type="SNP6_genomicSegment"
        if cancer =="COADREAD" and file in ["Gistic2_CopyNumber_Gistic2_all_data_by_genes",
                                            "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
                                            "Gistic2_CopyNumber_Gistic2_focal_data_by_genes",
                                            "mutation"]:            
            continue

        process (outDir, cancer, c1, c2, file, REALRUN,type)
    return


def process (outDir, cancer, c1, c2, file, REALRUN,type):
    if file =="":
        return

    c1dir = outDir+c1+"/"
    c2dir = outDir+c2+"/"

    c1file = c1dir+file
    c2file = c2dir+file

    if type in ["RPPA"]: #test file is the same
        os.system("cut -f 1 "+c1file +" >tmp1")
        os.system("cut -f 1 "+c2file +" >tmp2")
        SAME = filecmp.cmp("tmp1", "tmp2")

        if not SAME:
            print "files are not the same"
            return
        
        os.system("sort "+c1file +" >tmp1")
        os.system("sort "+c2file +" >tmp2")
        os.system("join -t$'\t' -1 1 -2 1 tmp1 tmp2 > tmp3")
        os.system("grep sample tmp3 > "+ outDir+cancer+"/"+file)
        os.system("grep -v sample tmp3 >> "+ outDir+cancer+"/"+file)
        
    elif os.path.exists(c1file) and os.path.exists(c2file):
        if REALRUN:
            if type !="SNP6_genomicSegment":
                os.system("cut -f 2- "+c2file +" >tmp")
                os.system("paste "+c1file+" tmp > "+outDir+cancer+"/"+file)
            else:
                os.system("cat "+c1file+" "+c2file +" > "+outDir+cancer+"/"+file)
                # segToMatrix
                outfile=outDir+cancer+"/"+file
                gMoutput = outDir+cancer+"/"+string.replace(file,"_genomicSegment","")
                os.system("python seg2matrix/segToMatrix.py "+outfile +" seg2matrix/refGene_hg18 "+ gMoutput)

    J={}
    iHandle =open(c1file+".json","r")
    Jinput= json.loads(iHandle.read())
    J[":dataSubType"]=  Jinput[":dataSubType"]
    J["version"]=  Jinput["version"]
    J["type"]=  Jinput["type"]
    J["cgDataVersion"]=     Jinput["cgDataVersion"]
    s= Jinput["shortTitle"]
    s = string.replace(s,c1,cancer)
    J["shortTitle"]= s
    J["label"] = J["shortTitle"] 
    J[":probeMap"]=Jinput[":probeMap"]
    J["anatomical_origin"]= TCGAUtil.anatomical_origin[cancer]
    J["sample_type"]=["tumor"]
    J["primary_disease"]=TCGAUtil.cancerGroupTitle[cancer]
    J["cohort"] ="TCGA_"+cancer
    J['domain']="TCGA"
    J['tags']=["cancer"]+ TCGAUtil.tags[cancer]
    J['owner']="TCGA"

    if Jinput.has_key("gdata_tags"):
        J["gdata_tags"]= Jinput["gdata_tags"]
    if Jinput.has_key("gain"):
        J["gain"]= Jinput["gain"]
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
        cfile = c1dir+gMoutput+".matrix.json"
        iHandle =open(cfile+".json","r")
        Jinput= json.loads(iHandle.read())
        J[":dataSubType"]=  Jinput[":dataSubType"]
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
