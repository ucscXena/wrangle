import string, sys, os
import csv
from itertools import izip
import uuid

def process (infile, outfile, input_delimiter, k):
    fout=open(outfile,'w')
    fout.close()
    fout=open(outfile,'a')

    #header
    fin=open(infile,'rU')
    line = fin.readline()
    totalN = len(string.split(line, input_delimiter))
    print (totalN)
    fin.close()

    count = 0
    tmpFile = str(uuid.uuid4())
    os.system("mkdir " + tmpFile)
    for i in range (0, totalN, k):
        count = count + 1
        tmpinfile = tmpFile + "/." + str(count)
        start = i + 1 #linux cut
        command = "cut -d '" + input_delimiter + "' -f " + str(start) + "-" + str(start+k-1) +" " + infile + " > " + tmpinfile
        os.system(command)

        fin=open(tmpinfile,'rU')
        a = izip(*csv.reader(fin, delimiter = input_delimiter))
        fin.close()
        csv.writer(fout, delimiter="\t", lineterminator="\n").writerows(a)
        os.system("rm -rf " + tmpinfile)

    fout.close()
    os.system("rm -rf " + tmpFile)

if len(sys.argv[:])!= 5:
    print ("python transpose.py fin fout input_delimiter(eg. , $\'\\t\') size(e.g.200 larger the file, smaller the size)")
    sys.exit()

infile = sys.argv[1]
outfile= sys.argv[2]
input_delimiter = sys.argv[3]
k = int(sys.argv[4])

process(infile, outfile, input_delimiter, k)

