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
    print ("python cluster.tsv.json.py dataDir")
    sys.exit()

dir = sys.argv[1]
configfile = dir + '/config'
metaDic = parse(configfile)

version = metaDic["version"]
cohort = metaDic["cohort"]

J ={}
J["type"] = "clinicalMatrix"
J["version"] = version
J["cohort"] = cohort
J["label"] = "louvain clusters"
J["dataSubtype"] = "phenotype"
J[":clinicalFeature"]="/HCA/clinicalFeature"
J["wrangling_procedure"]="Michael Krauss scanpy clustering run"

outfile = dir + '/cluster.tsv.json'
fout = open(outfile,'w')
json.dump(J, fout, indent =4)
fout.close()
