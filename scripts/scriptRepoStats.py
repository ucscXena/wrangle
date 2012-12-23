import sys,os,string,copy,json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")
from SampleMapNew import *
from CGDataLib import *
from ClinicalMatrixNew  import *
from ClinicalFeatureNew  import *
from CGDataUtil import *

inDir="data_flatten"
outfile ="statReport"

def genomicMatrixSampleNumber(file):
    fin = open(file,'r')
    n =len(string.split(fin.readline(),"\t"))-1
    fin.close()
    return n

def runRepoStats(inDir, outfile):
    dir = inDir
    bookDic={}
    ignore=1
    #repo
    bookDic=cgWalk(dir,ignore)
    if not bookDic :
        print "repo has problem"
        return 0
    #collcect stats
    statDic={}
    dataSubTypes=[]
    tDataset=0
    for dataset in bookDic:
        obj = bookDic[dataset]
        if obj["type"] !="genomicMatrix":
            continue
        name= obj["name"]
        groupTitle= obj["groupTitle"]
        dataSubType = obj[":dataSubType"]
        path=obj["path"]
        tDataset = tDataset +1
        if dataSubType not in dataSubTypes:
            dataSubTypes.append(dataSubType)
        if statDic.has_key(groupTitle):
            if statDic[groupTitle].has_key(dataSubType):
                statDic[groupTitle][dataSubType]= statDic[groupTitle][dataSubType]+genomicMatrixSampleNumber(path)
            else:
                statDic[groupTitle][dataSubType]= genomicMatrixSampleNumber(path)
        else:
            statDic[groupTitle]={}
            statDic[groupTitle][dataSubType]= genomicMatrixSampleNumber(path)

    #output stats
    fout=open(outfile,"w")
    dataSubTypes.sort()
    fout.write("data group\t"+string.join(dataSubTypes,"\t")+"\n")
    index = statDic.keys()
    index.sort()
    tSample=0
    for groupTitle in index:
        fout.write(groupTitle)
        for dataSubType in dataSubTypes:
            if statDic[groupTitle].has_key(dataSubType):
                fout.write("\t"+str(statDic[groupTitle][dataSubType]))
                tSample= tSample +statDic[groupTitle][dataSubType]
            else:
                fout.write("\t")
        fout.write("\n")
    fout.write("total dataset = "+ str(tDataset)+"\n")
    fout.write("total sample = "+ str(tSample)+"\n")
    fout.close()

runRepoStats(inDir, outfile)


