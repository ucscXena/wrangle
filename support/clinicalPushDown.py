import string, os, sys

def mapping (mfile):
    fin = open(mfile,'r')
    mapDic ={}
    for line in fin.readlines():
        data = line[:-1].split('\t')
        if len(data)==2:
            xenaSampleID, sourceID = data
            if sourceID not in mapDic:
                mapDic [sourceID]=[xenaSampleID]
            elif sourceID not in mapDic[sourceID]:
                mapDic [sourceID].append(xenaSampleID)

    fin.close()
    return mapDic

def pushDown (inputFile, outputFile, mapDic, _PATIENT_field_header):
    fin = open(inputFile, 'r')
    fout = open(outputFile, 'w')

    #header
    data = fin.readline().split('\t')
    N= len(data)
    fout.write("xena_sample\t" + _PATIENT_field_header + "\t")
    fout.write('\t'.join(data[1:]))

    #data
    for line in fin.readlines():
        data = line.split('\t')
        if (N!= len(data)):
            print (data, len(data))
        patientID = data[0]
        if patientID in mapDic:
            for sample in mapDic[patientID]:
                donor = data[0]
                fout.write(sample+'\t'+donor+'\t')
                fout.write('\t'.join(data[1:]))
                if (data[-1][-1]!='\n'):
                    fout.write('\n')

    fin.close()
    fout.close()

if len(sys.argv[:]) not in  [4, 5]:
    print ("python clinicalPushDown.py patientLevelFile mapping(sample patient no_header) sampleLevelFile _PATIENT_field_header(optional)")
    print
    sys.exit()

mappingfile = sys.argv[2]
mapDic = mapping (mappingfile)
inputfile = sys.argv[1]
outputfile = sys.argv[3]

if len(sys.argv[:]) ==5:
    _PATIENT_field_header = sys.argv[4]
else:
    _PATIENT_field_header = '_PATIENT'

pushDown (inputfile, outputfile, mapDic, _PATIENT_field_header)
