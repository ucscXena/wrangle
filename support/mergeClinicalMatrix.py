import sys,string,os
import json,datetime
import math
import copy
import uuid

def merge (infile_list, outfile):
    parameters=[]
    data ={}

    tmpList =[]
    for infile in infile_list:
        if not os.path.exists(infile):
            print infile, "not exists, skip"
        else:
            tmpList.append(infile)

    infile_list = tmpList

    for infile in infile_list:
        addData (infile, data, parameters)

    output (data, parameters, outfile)
    return

def addData (infile, DATA, parameters):
    fin=open(infile,'U')
    line = fin.readline()
    headers = string.split(line[:-1], '\t')

    for i in range (1, len(headers)):
        header = headers[i]
        if header not in parameters:
            parameters.append(header)
    
    for line in fin.readlines():
        items = string.split(line[:-1], '\t')
        sample = items[0]
        if sample not in DATA:
            DATA[sample] ={}
        for i in range (1, len(headers)):
            header = headers[i]
            if header in DATA[sample]:
                print "old:", DATA[sample][header], "new:", items[i]
            DATA[sample][header] = items[i]

    fin.close()
    return


def  output (data, parameters, outfile):
    fout=open(outfile,'w')

    fout.write('sample\t'+ string.join(parameters, '\t') + '\n')

    samples = data.keys()

    for sample in samples:
        fout.write(sample)
        for parameter in parameters:
            if parameter in data[sample]:
                fout.write('\t'+ data[sample][parameter])
            else:
                fout.write('\t')
        fout.write('\n')

    fout.close()
    return

if __name__ == "__main__" and len(sys.argv[:]) < 3:
    print "python mergeClinicalMatrix.py outfile infiles"
    sys.exit()

if __name__ == "__main__":
    outfile = sys.argv[1]
    infile_list =sys.argv[2:]
    print infile_list, outfile
    merge (infile_list, outfile)
