import sys,os

def copyOriginalFile (filename):
    os.system("cp " +filename +" "+filename+"_BK")
    return filename+"_BK"

def mapping (mfile):
    fin = open(mfile,'r')
    mapDic ={}
    for line in fin.readlines():
        data = line[:-1].split()
        xenaSampleID, sourceID = data[0:2]
        if sourceID not in mapDic:
            mapDic [sourceID]=xenaSampleID
        else:
            print ("ERROR", sourceID, "exists")
    fin.close()
    return mapDic

def switchIDs (inputFile, outputFile, mapDic, strict_inclusion = False, exclusion_mapDic = {}):
    fin = open(inputFile,'r')
    fout = open(outputFile,'w')

    #header
    fout.write(fin.readline())

    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = line[:-1].split('\t')
        id = data[0]
        if id not in mapDic:
            newID = id
            if strict_inclusion:
                print (id, "not found skip")
                continue
            else:
                print (id, "not found keep")
        else:
            newID = mapDic[id]

        if id in exclusion_mapDic:
            print (id, "skip")
            continue

        data[0]=newID
        fout.write('\t'.join(data)+'\n')
    fin.close()
    fout.close()

if len(sys.argv[:]) < 4 or  len(sys.argv[:]) > 5 :
    print
    print ("xenaIDreplaceMutationVector_clinicalMatrix.py mutationVector_clinMatrix_file(with_Header) \
        inclusion_mappingfile(new current) strict_inclusion(T/F) exclusion_mappping_file(optional)")
    print ("mapping file: no header")
    print ("              new_id old_id")
    print
    sys.exit()

inputfile = sys.argv[1]

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)

strict_inclusion = sys.argv[3]
if strict_inclusion not in ["T","F"]:
    print ("bad strict_inclusion value")
    sys.exit()
if strict_inclusion == "T":
    strict_inclusion = True
else:
    strict_inclusion = False

exclusion_mapDic = {}
if len(sys.argv[:]) == 5:
    exclusion_mappingfile = sys.argv[4]
    exclusion_mapDic = mapping (exclusion_mappingfile)

backupfile = copyOriginalFile(inputfile)
switchIDs (backupfile, inputfile, mapDic, strict_inclusion, exclusion_mapDic)
