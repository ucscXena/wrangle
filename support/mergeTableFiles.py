import sys,os
import uuid

def checkDataSetsExists(inputs):
    retList = []
    for file in inputs:
        if os.path.exists(file):
            retList.append(file)
        else:
            print file +" does not exist, skip"
    return retList

def mergeTableFiles (noDup, outfile, datasets):
    HEADER=""
    for path in datasets:
        fin = open(path,'r')
        header = fin.readline()
        if header[0]=="#":
            header = header[1:]
        if HEADER=="":
            HEADER= header
        if header != HEADER:
            print "error, different header\n"
            sys.exit()
        fin.close()

    if noDup == "noDup":
        tmpFile= str(uuid.uuid4())
        newTmpFile= str(uuid.uuid4())
        fdata = open(tmpFile,'w')
        for path in datasets:
            fin = open(path,'r')
            fin.readline()
            fdata.write(fin.read())
            fin.close()
        fdata.close()
        os.system("sort " + tmpFile +" | uniq > "+ newTmpFile )
        os.system("mv "+ newTmpFile + " " + tmpFile)

        fout=open(outfile,'w')
        fout.write(HEADER)
        fout.close()
        os.system("cat "+ tmpFile + " >> " + outfile)
        os.system("rm " + tmpFile)

    else:
        fout=open(outfile,'w')
        fout.write(HEADER)
        for path in datasets:
            fin = open(path,'r')
            fin.readline()
            fout.write(fin.read())
            fin.close()
        fout.close()


if len(sys.argv[:]) <4:
    print "python mergeTableFiles.py noDup/Dup outputfile inputfile(s)"
    sys.exit()

noDup = sys.argv[1]
if noDup not in ["noDup","Dup"]:
    print "python mergeTableFiles.py noDup/Dup outputfile inputfile(s)"
    sys.exit()

output= sys.argv[2]
inputs = sys.argv[3:]

inputs = checkDataSetsExists(inputs)
mergeTableFiles (noDup, output, inputs)

