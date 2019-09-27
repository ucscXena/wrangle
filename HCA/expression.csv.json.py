import sys, json

if len(sys.argv[:])!= 5:
    print "python expression.csv.json.py dataDir version cohort url"
    sys.exit()

dir = sys.argv[1]
version = sys.argv[2]
cohort = sys.argv[3]
url = sys.argv[4]

outfile = dir + '/expression.tsv.json'
fout = open(outfile,'w')

J ={}
J["type"] = "genomicMatrix"
J["version"] = version
J["cohort"] = cohort
J["label"] = "Optimus count"
J["unit"] = "count"
J["url"] = url
J["colNormalization"] = "log2(x)"
J["dataSubtype"] = "single cell RNAseq gene expression"

json.dump(J, fout, indent =4)
fout.close()