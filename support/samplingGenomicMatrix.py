import sys, os, string
import random, popen2

def dimention(infile):
    fin =open(infile,'r')
    line =fin.readline()
    column = len(string.split(line,'\t'))-1
    fin.close()

    r, w, e = popen2.popen3("wc -l "+ infile)
    row = int(string.split(r.read())[0])-1
    return row, column

def sampling(infile, nPerColumn):
    fin =open(infile,'r')
    line =fin.readline()
    retList=[]
    while 1:
        line = fin.readline()[:-1]
        if line =="":
            break
        values = string.split(line,'\t')[1:]
        l = random.sample(values, nPerColumn)
        for value in l:
            try:
                value = float(value)
                retList.append(value)
            except:
                continue
    fin.close()
    retList.sort()
    return retList

def output(ret_l, outfile):
    fout =open(outfile,'w')
    for value in ret_l:
        fout.write(str(value)+'\n')
    fout.close()

if len(sys.argv[:])!=4:
    print "python samplingGenomicMatrix.py matrixfile samplingN outfile"
    sys.exit()

infile = sys.argv[1]
N = int(sys.argv[2])
outfile =sys.argv[3]

row, column =dimention(infile)

if row *column <= N:
    print "use the whole file"
else:
    nPerColumn = N/row *2
    if nPerColumn ==0:
        nPerColumn =1
    sampleList = sampling(infile, nPerColumn)
    output(sampleList, outfile)
