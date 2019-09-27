import sys, json

if len(sys.argv[:])!= 4:
    print "python cells.csv.json.py dataDir version cohort"
    sys.exit()

dir = sys.argv[1]
version = sys.argv[2]
cohort = sys.argv[3]

outfile = dir + '/cells.tsv.json'
fout = open(outfile,'w')

J ={}
J["type"] = "clinicalMatrix"
J["version"] = version
J["cohort"] = cohort
J["label"] = "cell metadata"

json.dump(J, fout, indent =4)
fout.close()