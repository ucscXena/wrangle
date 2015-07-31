import sys,os,string,copy,json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")
from SampleMapNew import *
from CGDataLib import *
from ClinicalMatrixNew  import *
from ClinicalFeatureNew  import *
from CGDataUtil import *

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

        dataSubType = obj["dataSubType"]
        path=obj["path"]
        tDataset = tDataset +1
        if dataSubType not in dataSubTypes:
            dataSubTypes.append(dataSubType)
        if statDic.has_key(dataSubType):
            statDic[dataSubType]= statDic[dataSubType]+genomicMatrixSampleNumber(path)
        else:
            statDic[dataSubType]= genomicMatrixSampleNumber(path)

    #output stats
    dataSubTypes.sort()
    index = statDic.keys()
    index.sort()
    tSample=0

    for dataSubType in dataSubTypes:
        fout.write(dataSubType+"\t"+str(statDic[dataSubType]))
        tSample= tSample +statDic[dataSubType]
        fout.write("\n")

    fout.write("total dataset = "+ str(tDataset)+"\n")
    fout.write("total sample = "+ str(tSample)+"\n")


outfile ="statReport"
fout=open(outfile,"w")

fout.write("all data\n")
inDir="CAVM/"
runRepoStats(inDir, fout)

fout.write("\nopen data\n")
inDir="CAVM/public"
runRepoStats(inDir, fout)

fout.close()


