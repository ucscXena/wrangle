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
    print "python cells.csv.json.py dataDir"
    sys.exit()

dir = sys.argv[1]
configfile = dir + '/config'
metaDic = parse(configfile)

version = metaDic["version"]
cohort = metaDic["cohort"]
url = metaDic["url"]

outfile = dir + '/cells.tsv.json'
fout = open(outfile,'w')

J ={}
J["type"] = "clinicalMatrix"
J["version"] = version
J["cohort"] = cohort
J["label"] = "cell metadata"
J["dataSubtype"] = "phenotype"
J["url"] = url

json.dump(J, fout, indent =4)
fout.close()