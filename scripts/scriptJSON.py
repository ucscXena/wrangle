import os,sys, json
os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

from CGDataLib import *

def process(dir):
    bookDic=cgWalk(dir,1)
    for object in bookDic:
        J= bookDic[object]
        if J.has_key("tags"):
            if isinstance(J["tags"], basestring):
                tags = string.split(J["tags"],",")
                list =[]
                for tag in tags:
                    list.append(str(string.strip(tag)))
                J["tags"]=list
            else:
                list=[]
                for tag in J["tags"]:
                    list.append(str(string.strip(tag)))
                J["tags"]=list

        if J.has_key("sample_type"):
            if isinstance(J["sample_type"], basestring):
                sts = string.split(J["sample_type"],",")
                list =[]
                for st in sts:
                    list.append(str(string.strip(st)))
                J["sample_type"]=list
            else:
                list=[]
                for st in J["sample_type"]:
                    list.append(str(string.strip(st)))
                J["sample_type"]=list

        if J.has_key("sample_type"):
            if "cell line" in J["sample_type"]:
                if J.has_key("tags"):
                    if "cell lines" not in J["tags"]:
                        J["tags"].append("cell lines")
                else:
                    J["tags"]=["cell lines"]

        if J.has_key(":dataSubType"):
            if J[":dataSubType"]=="geneExp":
                if J.has_key("gdata_tags"):
                    if "transcription" not in J["gdata_tags"]:
                        J["gdata_tags"].append("transcription")
                else:
                    J["gdata_tags"]=["transcription"]

                if J.has_key("tags"):
                    if "transcription" in J["tags"]:
                        J["tags"].remove("transcription")

        if J.has_key(":dataSubType"):
            if J[":dataSubType"]=="geneExp":
                if J.has_key("gdata_tags"):
                    if "mRNA" not in J["gdata_tags"]:
                        J["gdata_tags"].append("mRNA")
                else:
                    J["gdata_tags"]=["mRNA"]

                if J.has_key("tags"):
                    if "mRNA" in J["tags"]:
                        J["tags"].remove("mRNA")

        if J["type"] in ["genomicMatrix","genomicSegment","sampleMap","clinicalMatrix"]:
            if J["type"] == "sampleMap":
                if J.has_key('cohort'):
                    J.pop('cohort')

            if J.has_key("tags"):
                if "cancer" not in J["tags"]:
                    J["tags"].append("cancer")
            else:
                J["tags"]=["cancer"]

        if J.has_key("anatomical_origin"):
            if isinstance(J["anatomical_origin"], basestring):
                origins = string.split(J["anatomical_origin"],",")
                list =[]
                for origin in origins:
                    list.append(str(string.strip(origin)))
                J["anatomical_origin"]=list
            if "brain" in J["anatomical_origin"] or \
               "Brain" in J["anatomical_origin"]:
                if J.has_key("tags"):
                    if "neural" not in J["tags"]:
                        J["tags"].append("neural")
                else:
                    J["tags"]=["neural"]

        if (J["type"] in ["genomicMatrix","genomicSegment","mutationVector"]):
            if not J.has_key("min") or not J.has_key("max"):
                if J.has_key("gain"):
                    J["min"]= -1.0/float(J["gain"])
                    J["max"]= 1.0/float(J["gain"])
                else:
                    J["min"]= -1.0
                    J["max"]= 1.0
        if J.has_key(":dataSubType"):
            if J[":dataSubType"]=="cna":
                J[":dataSubType"]="copy number"
            if J[":dataSubType"]=="DNAMethylation":
                J[":dataSubType"]="DNA methylation"
            if J[":dataSubType"]=="geneExp":
                J[":dataSubType"]="gene expression"
            if J[":dataSubType"]=="kinomeScreen":
                J[":dataSubType"]="kinase inhibition"
            if J[":dataSubType"]=="PARADIGM":
                J[":dataSubType"]="PARADIGM pathway activity"
            if J[":dataSubType"]=="protein":
                J[":dataSubType"]="protein expression RPPA"
            if J[":dataSubType"]=="somaticMutation":
                J[":dataSubType"]="somatic mutation"
            if J[":dataSubType"]=="siRNAViability":
                J[":dataSubType"]="cell viability"

        #temporary
        #if (J["type"] in ["genomicMatrix","genomicSegment","mutationVector"]) and J.has_key("shortTitle"):
        #    J["label"] = J["shortTitle"]

        fout=open(J['path']+".json","w")
        fout.write(json.dumps(J,indent=-1,sort_keys=True))
        fout.close()
            
dir="data_flatten/"
process(dir)
