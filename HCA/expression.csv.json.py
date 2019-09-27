import sys, json

def parse(configfile):
	fin = open(configfile, 'r')
	dic ={}
	for line in fin.readlines():
		key, value = line[:-1].split('\t')
		dic[key] = value
	fin.close()
	return dic

if len(sys.argv[:])!= 2:
    print ("python expression.csv.json.py dataDir")
    sys.exit()

dir = sys.argv[1]

configfile = dir + '/config'
metaDic = parse(configfile)

version = metaDic["version"]
cohort = metaDic["cohort"]
url = metaDic["url"]
exp_unit = metaDic["exp_unit"]
exp_label = metaDic["exp_label"]

J ={}
J["type"] = "genomicMatrix"
J["version"] = version
J["cohort"] = cohort
J["label"] = exp_label
J["unit"] = exp_unit
J["url"] = url
J["colNormalization"] = "log2(x)"
J["dataSubtype"] = "single cell RNAseq gene expression"

outfile = dir + '/expression.tsv.json'
fout = open(outfile,'w')
json.dump(J, fout, indent =4)
fout.close()

outfile = dir + '/expression.xena.tsv.json'
fout = open(outfile,'w')
json.dump(J, fout, indent =4)
fout.close()