import sys, os

def listing (mfile):
    fin = open(mfile,'r')
    list =[]
    for line in fin.readlines():
        # data = string.strip(string.split(line,'\t')[0])
        data = line.split('\t')[0].strip()
        if data not in list:
            list.append(data)
    fin.close()
    return list

def keepIDs (inputFile, outputFile, keep_list):
    goodPos=[]

    fin = open(inputFile,'r')
    ids = fin.readline()[:-1].split('\t')
    fin.close()

    for i in range (len(keep_list)):
        keepID = keep_list[i]
        try:
            index = ids.index(keepID)
            goodPos.append(index)
        except:
            continue

    fin = open(inputFile,'r')
    fout = open(outputFile, 'w')
    while 1:
        line = fin.readline()
        if line == '':
            break
        data = line[:-1].split('\t')
        good_data =[]
        #row name
        fout.write(data[0]+'\t')
        # rest of the line
        for pos in goodPos:
            good_data.append(data[pos])
        fout.write('\t'.join(good_data)+'\n')
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print ("python keepIDGenomicMatrix.py matrixInfile matrixOutfile(same order as in keep_list) keep_list(one_sample_id_per_line, first column ID, no header)")
    sys.exit(1)

listfile = sys.argv[3]
keep_list  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

keepIDs (inputfile, outputfile, keep_list)
