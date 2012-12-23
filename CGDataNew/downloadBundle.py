import sys,os,string,copy,json,re,datetime
from SampleMapNew import *
from CGDataLib import *
from ClinicalMatrixNew  import *
from ClinicalFeatureNew  import *
from CGDataUtil import *
ASSEMBLY ="hg18"

tarballDir = "TestdownloadTarball/"
outDir = "Testdownload/"
ROOT ="data_flatten/"

def probeMapWalk(probeMap):
    # all probeMap
    probeDic=cgWalk(probeMap,0)
    for dataset in probeDic:
        type = probeDic[dataset]['type']
        path = probeDic[dataset]['path']
        version = probeDic[dataset]['version']
        assembly = probeDic[dataset][":assembly"]
        
        if probeDic[dataset].has_key("group") and  probeDic[dataset]["group"] is not None:
            group = probeDic[dataset]["group"]
    return probeDic

def findProbe(probeDic, name, assembly):
    if name is None or name =="":
        return None
    if probeDic.has_key(name):
        if probeDic[name][":assembly"]== assembly:
            return probeDic[name]

    for dataset in probeDic:
        if probeDic[dataset].has_key("group") and \
           probeDic[dataset]["group"] ==name:
            if probeDic[dataset][":assembly"]== assembly:
                return probeDic[dataset]
    return None
    
def downloadBundle(inDir, PROBEMAP, testVersionOnly):
    dir = inDir
    #test inDir must be from data_flatten
    #if dir[0:12]!="data_flatten":
    #    print "wrong input data dir"
    #    return 1
    secureDir =string.replace(inDir, ROOT,"")
    secureDir =string.lower(string.split(secureDir,'/')[0])
    print secureDir
    #skip data_flatten/common/ and data_flatten/probeMap
    if secureDir in ["probemap","common"]:
        return 0
    #dir related to security    
    if secureDir not in ["gray","ispy","pancreas","public","su2c"]:  #,"tcga"]:
        print "wrong security dir", secureDir
        return 1
    
    bookDic={}
    ignore=0
    bookDic=cgWalk(dir,ignore)
    if not bookDic :
        print "repo has problem"
        return 1

    probeMapDic={}

    if not testVersionOnly:
        probeMapDic=probeMapWalk(PROBEMAP)

    badVersion=0
    for dataset in bookDic:
        type = bookDic[dataset]['type']
        if type != "genomicMatrix":
            continue

        if bookDic[dataset].has_key("redistribution") and bookDic[dataset]["redistribution"] ==False:
            print "skip redistribution=false"
            continue
        path = bookDic[dataset]['path']

        #version
        for clinDataset in bookDic:
            if  bookDic[clinDataset]['type'] !="clinicalMatrix" or bookDic[clinDataset][':sampleMap'] !=bookDic[dataset][':sampleMap']:
                continue
            cVersion = bookDic[clinDataset]['version']
            break
        gVersion = bookDic[dataset]['version']

        dG= makeDate(gVersion)
        dC= makeDate(cVersion)
        if dG<dC:
            version = cVersion
        else:
            version= gVersion
            
        #test version
        if len(version)!=10:
            print "bad version", dataset, version
            badVersion=1
            continue
        pattern= re.compile("[0-9]{4}-[0-9]{2}-[0-9]{2}")
        r = pattern.match(version)
        if r ==None:
            print "bad version", dataset, version
            badVersion=1
            continue

        if testVersionOnly:
            continue
    
        subdir = secureDir+"/"+dataset+"-"+version+"/"
        targetdir = dataset+"-"+version+"/"
        if os.path.exists( outDir+ subdir ):
            os.system("rm -rf "+outDir+subdir+"*")
        else:
            os.system("mkdir "+outDir+subdir)
        #os.system("rm -rf "+tarballDir + dataset +"-[0-9]{4}-[0-9]{2}-[0-9]{2}.tar.gz")

        print "dataset=", dataset

        #probeMap
        if type =="genomicMatrix" and not bookDic[dataset].has_key(":genomicSegment"):
            name = bookDic[dataset][":probeMap"]
            probeMap = findProbe(probeMapDic, name, ASSEMBLY)
            if probeMap == None:
                print "skip\n"
                continue
            path = probeMap['path']
            print "probeMap=", name, path
            os.system("cp "+ path +" "+outDir+ subdir)
            os.system("cp "+ path +".json "+outDir + subdir)

        #genomic data
        if bookDic[dataset].has_key(":genomicSegment"):
            dataset = bookDic[dataset][":genomicSegment"]
            path = bookDic[dataset]['path']
            type = bookDic[dataset]['type']
        else:
            path = bookDic[dataset]['path']
        os.system("cp "+ path +" "+outDir+subdir)
        os.system("cp "+ path +".json "+outDir+subdir)
        
        #sampleMap
        sampleMap = bookDic[dataset][':sampleMap']
        path = bookDic[sampleMap]['path']
        os.system("cp "+ path +" "+outDir+ subdir)
        os.system("cp "+ path +".json "+outDir + subdir)


        #clinicalMatrix
        for clinDataset in bookDic:
            if  bookDic[clinDataset]['type'] !="clinicalMatrix" or bookDic[clinDataset][':sampleMap'] !=sampleMap:
                continue
            print "clinicalMatrix=", bookDic[clinDataset]['name']
            path=bookDic[clinDataset]['path']
            os.system("cp "+ path +" "+outDir+ subdir)
            os.system("cp "+ path +".json "+outDir + subdir)

            #clinicalFeature
            if bookDic[clinDataset].has_key(":clinicalFeature"):
                print "clinicalFeature=", bookDic[clinDataset][":clinicalFeature"]
                name = bookDic[clinDataset][":clinicalFeature"]
                path = bookDic[name]['path']
                os.system("cp "+ path +" "+outDir+ subdir)
                os.system("cp "+ path +".json "+outDir + subdir)

        #md5
        currdir = os.path.abspath(os.path.curdir)
        os.chdir(outDir+subdir)
        os.system("md5sum * > md5.txt")
        os.chdir(currdir)
        
        #tar.gz
        if not os.path.exists(tarballDir +secureDir):
            os.system("mkdir "+tarballDir +secureDir)
        os.system("tar -czPf "+ tarballDir +"/"+subdir[:-1]+".tgz -C "+ outDir+secureDir+" "+targetdir)

        if badVersion:
            print "bad version number detected"

        print

    return badVersion    
