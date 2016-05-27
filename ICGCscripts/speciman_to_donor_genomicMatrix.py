import string, sys,os, uuid

def copyOriginalFile (filename):
    os.system("cp " +filename +" "+filename+"_BK")
    return filename+"_BK"

#old new
#old MUST be uniq
#new not uniq
def reverse_mapping (mfile): 
    fin = open(mfile,'U')
    mapDic_oldTOnew ={}
    mapDic_newToOld={}
    mapDic_oldSampleType ={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        oldId, newId, sample_type = data[0:3]
        if oldId  not in mapDic_oldTOnew:
            mapDic_oldTOnew [oldId]= newId 
            mapDic_oldSampleType[oldId] = sample_type
        else:
            print "ERROR", oldId, "exists"
        if newId not in mapDic_newToOld:
            mapDic_newToOld[newId]=[oldId]
        else:
            mapDic_newToOld[newId].append(oldId)
    fin.close()
    return mapDic_oldTOnew, mapDic_newToOld, mapDic_oldSampleType


#old new
#no normals
#if still multiple tumors, only pick primary tumors
def processBySampleType (current_samples, mapDic_oldTOnew, mapDic_newToOld, mapDic_oldSampleType):
    isNormal = lambda x : string.strip(string.split(x,"-")[0]) in ["Normal","normal"]

    dic ={} 
    keepPos =[0]
    keepPosWithDup=[]
    dupDic ={} #newid: [pos,pos,..]

    for i in range(1,len(current_samples)):
        sample = current_samples[i]
        if isNormal(mapDic_oldSampleType[sample]):
            continue
        newid = mapDic_oldTOnew[sample]
        if newid not in dic:
            dic[newid]=[i]
        else:
            dic[newid].append(i)

    for id in dic.keys():
        if len(dic[id])!=1:
            tmp_list=[]
            for i in dic[id]:
                type = mapDic_oldSampleType[ samples[i]]
                if string.find(type,"Primary tumour")!=-1:
                    tmp_list.append(i)
            dic[id]=tmp_list
        if len(dic[id])==1:
            keepPos.extend(dic[id])
        else:
            keepPosWithDup.extend(dic[id])
            dupDic[id]=dic[id]
        #print id, dic[id]
        #for i in dic[id]:
        #    print samples[i], mapDic_oldSampleType[ samples[i]]

    keepPos.sort()
    keepPosWithDup.sort()
    print len(keepPos), len(keepPosWithDup), len(current_samples)
    return keepPos, keepPosWithDup, dupDic
    

def unixSimple(inputfile, keepPos):
    tmpfile = str(uuid.uuid4()) 
    s= "cut -f "
    for i in keepPos[:-1]:
        s= s+str(1+i)+","
    s = s + str(1+keepPos[-1]) + " " +inputfile +" > " + tmpfile
    os.system(s)
    return tmpfile

def samplesInFile (inputFile):
    fin = open(inputFile,'U')
    ids = string.split(fin.readline()[:-1],'\t')
    fin.close()
    return ids


def switchFirstlineIDs (inputFile, outputFile, mapDic):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    
    ids = string.split(fin.readline()[:-1],'\t')
    newIds =["donor"]
    for id in ids[1:]:
        if not id in mapDic:
            print id
            xenaID = id
        else:
            xenaID = mapDic[id]
        newIds.append(xenaID)
    fout.write(string.join(newIds,'\t')+'\n')
    while 1:
        line = fin.readline()
        if line =="":
            break
        fout.write(line)
    fin.close()
    fout.close()


if len(sys.argv[:])!=4:
    print 
    print "python speciman_to_donor_genomicMatrix.py matrixfile mappingfile(current new sampleType)"
    print 
    print "mapping file: no header"
    print "              current_id new_id sample_type"
    print 
    print "This script only generate these donors:"
    print "1. simple tumor (could be any type of tumor)"
    print "2. if a donor has multiple tumor, but with only one primary tumor, we use the primary tumor data"
    print "the rest of donor, we ignore."
    print 

    sys.exit()

mappingfile = sys.argv[2]
inputfile = sys.argv[1]
samples = samplesInFile(inputfile)
mapDic_oldTOnew, mapDic_newToOld, mapDic_oldSampleType = reverse_mapping (mappingfile)
keepPos, keepPosWithDup, dupDic  = processBySampleType (samples, mapDic_oldTOnew, mapDic_newToOld, mapDic_oldSampleType)

tmpfile = unixSimple(inputfile, keepPos)
outputfile = sys.argv[3]
switchFirstlineIDs (tmpfile, outputfile, mapDic_oldTOnew)

print tmpfile

