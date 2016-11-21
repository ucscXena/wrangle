import json,re,subprocess
import os,sys,copy,string,datetime
from ClinicalMatrixNew  import *
from SampleMapNew import *
from CGDataUtil import *

def goodLetterName(state):
    regex = re.compile('[a-z,A-Z,0-9,"_"]')
    for i in range (0,len(state)):
        if not regex.search(state[i]):
            return False
    return True

def cgWalk(dir, ignore):
    bookDic={}
    for root,dir,files in os.walk(dir,followlinks=True):
        if root[-1]=="/":
            root =root[:-1]
        if files ==[]:
            continue
        for file in files:
            if file[-5:]!= ".json":
                continue
            base =file [:-5]
            if os.path.getsize(root+"/"+file) ==0:
                if ignore:
                    continue
                else:
                    print "json file size 0", root+"/"+file
                    return 0

            #load .json
            J = json.loads(open(root+"/"+base+".json").read())
            J['path']= root+"/"+base
            if J.has_key("name"):
                name = J['name']
            else:
                name = root+"/"+base
                J["name"]= name
            if not J.has_key('type'):
                J['type']=None

            if J['type']!=None and not os.path.exists(root+"/"+base):
                if ignore:
                    continue
                else:
                    if J['type']=="cohort":
                        continue
                    print "data file does not exist", root+"/"+base
                    continue
                    #return 0

            #json checking if valid
            if bookDic.has_key(name):
                if ignore:
                    pass
                else:
                    print "CRITICAL ERROR already has this name in repo", name
                    return 0
            if type in ["genomicMatrix","clinicalMatrix","mutationVector"] and not trackName_good(name): 
                print "ERROR name has bad characters or too long", name
                return 0
            bookDic[name]=J
    return bookDic
                
def collectNamesBelongToSampleMap(bookDic, sampleMapName):
    list=[]
    for name in bookDic.keys():
        J = bookDic[name]
        if J.has_key(':sampleMap') and J[':sampleMap'] == sampleMapName:
            list.append(name)
    return list

def collectProbeMapBelongToSampleMap(bookDic, sampleMapName):
    probeMaps= []
    for dataset in collectNamesBelongToSampleMap(bookDic, sampleMapName):
        if bookDic[dataset].has_key(":probeMap"):
            probeMaps.append(bookDic[dataset][":probeMap"])
    list=[]
    for name in bookDic.keys():
        J = bookDic[name]
        if J["type"]=="probeMap" and name in probeMaps:
            list.append(name)
    return list

def collectSampleMaps(bookDic):
    sampleMaps={}
    for name in bookDic.keys():
        J = bookDic[name]
        if J['type'] == "sampleMap":
            if name not in sampleMaps:
                sampleMaps[name]=[]
    for name in bookDic.keys():
        J = bookDic[name]
        if J.has_key(':sampleMap') and sampleMaps.has_key(J[':sampleMap']):
            sampleMaps[J[':sampleMap']].append(name)
    return sampleMaps

def collectMissingSampleMaps(bookDic):
    sampleMaps={}
    missingMaps={}
    for name in bookDic.keys():
        J = bookDic[name]
        if J['type'] == "sampleMap":
            sName =J['name']
            if sName not in sampleMaps:
                sampleMaps[sName]=[]

    for name in bookDic.keys():
        J = bookDic[name]
        if J.has_key(':sampleMap') and not sampleMaps.has_key(J[':sampleMap']):
            if J[':sampleMap'] not in missingMaps:
                missingMaps[J[':sampleMap']]=[]
    for name in bookDic.keys():
        J = bookDic[name]
        if J.has_key(':sampleMap') and missingMaps.has_key(J[':sampleMap']):
            missingMaps[J[':sampleMap']].append(name)
    return missingMaps


def cgDataMergeJSON(J1, J2, name):
    J3 = {}
    for key in J1.keys():
        if key =="url":
            J3["url"]= J1["url"]
        if key =="cohort":
            J3["cohort"]= J1["cohort"]

    for key in J2.keys():
        if key =="cohort":
            J3["cohort"]= J2["cohort"]
        if key =="url":
            if J3.has_key("url") and string.find(J3["url"],J2["url"])==-1:
                J3["url"]= J3["url"]+", " +J2["url"]
            else:
                J3["url"]= J2["url"]

    J3['name']=name
    J3['type']="clinicalMatrix"
    J3[':sampleMap']=J1[":sampleMap"]
    J3["dataSubType"] = "phenotype"
    return J3

def checkIdsAllIn(sMap, bookDic):
    sMapChanged=0
    for name in bookDic.keys():
        obj = bookDic[name]
        if obj.has_key(':sampleMap') and obj[':sampleMap']==sMap.getName():
            if obj['type'] in ["clinicalMatrix","mutationVector"]:
                #get matrix obj
                path = obj['path']
                fin=open(path,'r')
                fin.readline()
                for line in fin.readlines():
                    if string.strip(line)=="":
                        break
                    sample = string.split(line,'\t')[0]
                    if not sMap.inData(sample):
                        sMap.addNode(sample)
                        sMapChanged=1
                        
            elif obj['type']=="genomicMatrix":
                path = obj['path']
                fin = open(path,'r')
                for sample in string.split(string.strip(fin.readline()),'\t')[1:]:
                    if not sMap.inData(sample):
                        sMap.addNode(sample)
                        sMapChanged=1
                fin.close()
            else:
                print obj['type'],"need to write code for this"
    return sMapChanged


def getAllGenomicIds(sMap, bookDic):
    allSamples=[]
    for name in bookDic.keys():
        obj = bookDic[name]
        if obj.has_key(':sampleMap') and obj[':sampleMap']==sMap.getName():
            if obj['type']=="clinicalMatrix":
                continue
            elif obj['type']=="genomicMatrix":
                path = obj['path']
                fin = open(path,'r')
                for sample in string.split(string.strip(fin.readline()[:-1]),'\t')[1:]:
                    if sample not in allSamples:
                        allSamples.append(sample)
                fin.close()
            elif obj['type']=="mutationVector":
                path = obj['path']
                fin=open(path,'r')
                fin.readline()
                for line in fin.readlines():
                    if string.strip(line)=="":
                        break
                    sample = string.split(line,'\t')[0]
                    if sample not in allSamples:
                        allSamples.append(sample)
            else:
                print obj['type'],"need to write code for this"
    return allSamples

def makeDate(isoFormatDate):
    data = string.split(isoFormatDate,"-")
    if len(data)!=3:
        return None
    for i in range (0,len(data)):
        try:
            int(data[i])
        except TypeError:
            return None
    y,m,d = data
    DATE = datetime.date(int(y), int(m), int(d))
    return DATE


