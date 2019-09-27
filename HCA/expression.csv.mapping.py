import sys, os

def mapping(infile, outfile):
	os.system("cut -d , -f 1 " + infile + ' | tail -n + 2 | cat -n > ' + outfile)

if len(sys.argv[:])!= 2:
    print "python expression.csv.mapping.py dataDir"
    sys.exit()

dir = sys.argv[1]

infile = dir + '/expression.csv'
outfile= dir + '/mapping'

mapping(infile, outfile) 