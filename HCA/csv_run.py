import sys, os

def parse(configfile):
	fin = open(configfile, 'r')
	dic ={}
	for line in fin.readlines():
		key, value = line[:-1].split('\t')
		dic[key] = value
	fin.close()
	return dic

if len(sys.argv[:])!= 3:
	print ("python csv_run.py dataDir runExpression(0,1)")
	sys.exit()
    
dir = sys.argv[1]
expressionRun = int(sys.argv[2])

configfile = dir + '/config'
if not os.path.exists(configfile):
	print ("not config file")
	sys.exit()

metaDic = parse(configfile)
if "version" in metaDic:
	version = metaDic["version"]
else:
	print ("missing version")
	sys.exit()

if "cohort" in metaDic:
	cohort = metaDic["cohort"]
else:
	print ("missing cohort")
	sys.exit()

if "url" in metaDic:
	url = metaDic["url"]
else:
	print ("missing url")
	sys.exit()

if "size" in metaDic:
	size = metaDic["size"]
else:
	print ("missing size")
	sys.exit()

if "exp_unit" in metaDic:
	size = metaDic["exp_unit"]
else:
	print ("missing exp_unit")
	sys.exit()

if "exp_label" in metaDic:
	size = metaDic["exp_label"]
else:
	print ("missing exp_label")
	sys.exit()

codedir = os.path.dirname(sys.argv[0])
os.system("python " + codedir + "/cells.csv.json.py " + dir)
os.system("python " + codedir + "/expression.csv.json.py " + dir)
os.system("python " + codedir + "/cells.csv.py " + dir)
if (expressionRun):
	os.system("python " + codedir + "/expression.csv.py " + dir + " " + size)
