import string, sys,os

def copyOriginalFile (filename):
    os.system("cp " +filename +" "+filename+"_BK")
    return filename+"_BK"

#new old
#old MUST be uniq
def mapping (mfile):  
    fin = open(mfile,'U')
    mapDic ={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        xenaSampleID, sourceID = data[0:2]
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
    print 
    print "python xenaIDreplaceGenomicMatrix.py matrixfile mappingfile(new current)"
    print "mapping file: no header"
    print "              new_id old_id"
    print 
    sys.exit()

mappingfile = sys.argv[2]
inputfile = sys.argv[1]

backupfile = copyOriginalFile(inputfile)

mapDic = mapping (mappingfile)
switchFirstlineIDs (backupfile, inputfile, mapDic)

