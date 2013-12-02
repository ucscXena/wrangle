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

        if J["type"] in ["genomicMatrix","genomicSegment","sampleMap"]:
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

        #temporary
        if (J["type"] in ["genomicMatrix","genomicSegment"]) and J.has_key("shortTitle"):
            J["label"] = J["shortTitle"]

        fout=open(J['path']+".json","w")
        fout.write(json.dumps(J,indent=-1,sort_keys=True))
        fout.close()
            
dir="data_flatten/public/"
process(dir)
