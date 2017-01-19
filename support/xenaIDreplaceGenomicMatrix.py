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

def switchFirstlineIDs (inputFile, outputFile, mapDic, strict_inclusion = False, exclusion_mapDic = {}):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    ids = string.split(fin.readline()[:-1],'\t')
    newIds =[ids[0]]
    exclusion_pos = []
    inclusion_pos = [] #strict situation

    for i in range (1, len(ids)):
        id = ids[i]
        if not id in mapDic:
            print id, "not found"
            xenaID = id
        else:
            xenaID = mapDic[id]
            if strict_inclusion:
                inclusion_pos.append(i)

            if xenaID in newIds:
                print id, xenaID, "duplicate"

        if id in exclusion_mapDic:
            print id, "skip"
            exclusion_pos.append(i)
            
        newIds.append(xenaID)

    if len(exclusion_pos) == 0 and len(inclusion_pos) == 0:
        fout.write(string.join(newIds,'\t')+'\n')
    else:
        fout.write(newIds[0])
        for i in range(1, len(newIds)):
            if i not in exclusion_pos:
                if not strict_inclusion or  i in inclusion_pos:
                    fout.write('\t'+newIds[i])
        fout.write('\n')

    while 1:
        line = fin.readline()
        if line =="":
            break
        if len(exclusion_pos) == 0 and len(inclusion_pos) == 0:
            fout.write(line)
        else:
            data = string.split(line[:-1],'\t')
            fout.write(data[0])
            for i in range(1, len(data)):
                if i not in exclusion_pos:
                    if not strict_inclusion or  i in inclusion_pos:
                        fout.write('\t'+data[i])
            fout.write('\n')

    fin.close()
    fout.close()

if len(sys.argv[:]) < 4 or  len(sys.argv[:]) > 5 :
    print 
    print "python xenaIDreplaceGenomicMatrix.py matrixfile inclusion_mappingfile(new current) strict_inclusion(T/F) exclusion_mappping_file(optional)"
    print "mapping file: no header"
    print "              new_id old_id"
    print 
    sys.exit()

inputfile = sys.argv[1]

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)

strict_inclusion = sys.argv[3]
if strict_inclusion not in ["T","F"]:
    print "bad strict_inclusion value"
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

switchFirstlineIDs (backupfile, inputfile, mapDic, strict_inclusion, exclusion_mapDic)

