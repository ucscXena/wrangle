import string, sys,os

def copyOriginalFile (filename):
    os.system("cp " +filename +" "+filename+"_BK")
    return filename+"_BK"

def mapping (mfile):
    fin = open(mfile,'U')
    mapDic ={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        if len(data)==2:
            xenaSampleID, sourceID = data
            if sourceID not in mapDic:
                mapDic [sourceID]=xenaSampleID
            else:
                print "ERROR", sourceID, "exists"
    fin.close()
    return mapDic

def switchFirstlineIDs (inputFile, outputFile, mapDic):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    
    ids = string.split(fin.readline()[:-1],'\t')
    newIds =[ids[0]]
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

if len(sys.argv[:])!=3:
    print "xenaIDreplaceGenomicMatrix.py matrixfile mappingfile(new current)"
    sys.exit()

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)
inputfile = sys.argv[1]
backupfile = copyOriginalFile(inputfile)
switchFirstlineIDs (backupfile, inputfile, mapDic)
