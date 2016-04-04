import string, sys,os

def listing (mfile):
    fin = open(mfile,'r')
    list =[]
    for line in fin.readlines():
        data = string.strip(line)
        if data not in list:
            list.append(data)
    fin.close()
    return list

def keepIDs (inputFile, outputFile, keep_list):
    fin = open(inputFile,'U')
    
    goodPos=[]
    ids = string.split(fin.readline()[:-1],'\t')
    for i in range (0, len(ids)):
        if ids[i] in keep_list:
            goodPos.append(i)
            print ids[i], i
    fin.close()

    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line[:-1],'\t')
        fout.write(data[0])
        for i in range (1, len(data)):
            if i in goodPos:
                fout.write("\t"+data[i])
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python keepIDGenomicMatrix_slow.py matrixfile output remove_list(one_id_per_line)"
    sys.exit()

listfile = sys.argv[3]
keep_list  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

keepIDs (inputfile, outputfile, keep_list)
