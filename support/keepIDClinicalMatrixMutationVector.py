import string, sys,os

def listing (mfile):
    fin = open(mfile,'r')
    list =[]
    for line in fin.readlines():
        data = string.strip(string.split(line,'\t')[0])
        if data not in list:
            list.append(data)
    fin.close()
    return list

def keepIDs (inputFile, outputFile, keep_list):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    fout.write(fin.readline())
    for line in fin.readlines():
        data = string.split(line,'\t')
        id = data[0]
        if id in keep_list:
            fout.write(line)
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "keepIDClinicalMatrix.py matrixfile output keep_list(first_column_id)"
    sys.exit()

listfile = sys.argv[3]
keep_list  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

keepIDs (inputfile, outputfile, keep_list)
