import csv, sys
from itertools import izip

if len(sys.argv[:])!= 3:
    print "python transpose.py fin fout"
    sys.exit()
    
infile = sys.argv[1]
fin=open(infile,'rU')
a = izip(*csv.reader(fin,delimiter="\t"))

outfile= sys.argv[2]
fout=open(outfile,'w')
csv.writer(fout,delimiter="\t").writerows(a)

fout.close()
fin.close()
                
