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

def switchIDs (inputFile, outputFile, mapDic):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    #header
    fout.write(fin.readline())
    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line,'\t')
        id = data[0]
        if  id not in mapDic:
            newID = id
            print id
        else:
            newID = mapDic[id]
        data[0]=newID
        fout.write(string.join(data,'\t'))
    fin.close()
    fout.close()

if len(sys.argv[:])!=3:
    print "xenaIDreplaceMutationVector.py mutationVectorfile(with_Header) mappingfile(new current)"
    sys.exit()

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)
inputfile = sys.argv[1]
backupfile = copyOriginalFile(inputfile)
switchIDs (backupfile, inputfile, mapDic)
