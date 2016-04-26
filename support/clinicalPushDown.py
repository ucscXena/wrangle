import string,os,sys

def mapping (mfile):
    fin = open(mfile,'U')
    mapDic ={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        if len(data)==2:
            xenaSampleID, sourceID = data
            if sourceID not in mapDic:
                mapDic [sourceID]=[xenaSampleID]
            elif sourceID not in mapDic[sourceID]:
                mapDic [sourceID].append(xenaSampleID)

    fin.close()
    return mapDic

def pushDown (inputFile, outputFile, mapDic):
    fin = open(inputFile, 'U')
    fout = open(outputFile, 'w')

    #header
    data = string.split(fin.readline(),'\t')
    N= len(data)
    fout.write("xena_sample\t_PATIENT\t")
    fout.write(string.join(data[1:],'\t'))

    #data
    for line in fin.readlines():
        data = string.split(line,'\t')
        if (N!= len(data)):
            print data, len(data)
        patientID = data[0]
        if patientID in mapDic:
            for sample in mapDic[patientID]:
                donor = data[0]
                fout.write(sample+'\t'+donor+'\t')
                fout.write(string.join(data[1:],'\t'))
        else:
            fout.write(line)

    fin.close()
    fout.close()

if len(sys.argv[:])!= 4:
    print "python clinicalPushDown.py patientLevelFile mapping(sample patient) sampleLevelFile"
    sys.exit()

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)
inputfile = sys.argv[1]
outputfile = sys.argv[3]
pushDown (inputfile, outputfile, mapDic)
