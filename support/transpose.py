import string, sys, os
import csv
from itertools import izip
import uuid

def process (infile, outfile):
    fout=open(outfile,'w')
    fout.close()
    fout=open(outfile,'a')

    #header
    fin=open(infile,'rU')
    line = fin.readline()
    totalN = len(string.split(line,'\t'))
    fin.close()

    k = 200
    count = 0
    tmpFile = str(uuid.uuid4())
    for i in range (0, totalN, k):
        count = count + 1
        tmpinfile = tmpFile + "." + str(count)
        start = i + 1 #linux cut
        command = "cut -f " + str(start) + "-" + str(start+k-1) +" " + infile + " > " + tmpinfile
        os.system(command)

        fin=open(tmpinfile,'rU')
        a = izip(*csv.reader(fin,delimiter="\t"))
        fin.close()
        csv.writer(fout,delimiter="\t").writerows(a)
    fout.close()
    os.system("rm " + tmpFile + ".*")

if len(sys.argv[:])!= 3:
    print "python transpose.py fin fout"
    sys.exit()

infile = sys.argv[1]
outfile= sys.argv[2]

process(infile, outfile)

