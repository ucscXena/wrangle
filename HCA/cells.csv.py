import sys

if len(sys.argv[:])!= 2:
    print "python cells.csv.py dataDir"
    sys.exit()

dir = sys.argv[1]
infile = dir + '/cells.csv' 
outfile = dir + '/cells.tsv'
os.system("cat " + infile + " | tr , '\t' > " + ouptut)
