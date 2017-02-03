import string, sys
import uuid

import csv
from itertools import izip


def process (infile, outfile):
    fout=open(outfile,'w')

    tmpDir = str(uuid.uuid4())

    #header
    fin=open(infile,'rU')
    line = fin.readline()
    totalN = len(string.split(line,'\t'))
    fin.close()

    N =50
    count = 0
    for i in range (0, totalN, N):
        count = count + 1
        start = i + 1 #linux cut
        end = i + N #linux cut
        tmpoutput = "file_" + str(count)
        os.system("cut -f " + str(start) + "-" + str(end) + " " + infile + " > " + tmpDir + "/" + tmpoutput)

        ftmp=open(tmpoutput,'rU')
        a = izip(*csv.reader(ftmp,delimiter="\t"))
        csv.writer(fout,delimiter="\t").writerows(a)
        ftmp.close()

    fout.close()

if len(sys.argv[:])!= 3:
    print "python transpose_memEfficient.py fin fout"
    sys.exit()

infile = sys.argv[1]
outfile= sys.argv[2]

process(infile, outfile)
