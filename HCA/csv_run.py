import sys, os

def parse(configfile):
	fin = open(configfile, 'r')
	dic ={}
	for line in fin.readlines():
		key, value = line[:-1].split('\t')
		dic[key] = value
	fin.close()
	return dic

if len(sys.argv[:])!= 4:
	print ("python csv_run.py dataDir mappingRun(0,1) runExpression(0,1)")
	sys.exit()
    
dir = sys.argv[1]
mappingRun = int(sys.argv[2])
expressionRun = int(sys.argv[3])

configfile = dir + '/config'
if not os.path.exists(configfile):
	print ("not config file")
	sys.exit()

metaDic = parse(configfile)
if "version" not in metaDic:
	print ("missing version")
	sys.exit()

if "cohort" not in metaDic:
	print ("missing cohort")
	sys.exit()

if "url" not in metaDic:
	print ("missing url")
	sys.exit()

if "exp_unit" not in metaDic:
	print ("missing exp_unit")
	sys.exit()

if "exp_label" not in metaDic:
	print ("missing exp_label")
	sys.exit()

if "size" in metaDic:
	size = metaDic["size"]
else:
	print ("missing size")
	sys.exit()

codedir = os.path.dirname(sys.argv[0])
if mappingRun:
        os.system("python " + codedir + "/expression.csv.mapping.py " + dir)
os.system("python " + codedir + "/cells.csv.json.py " + dir)
os.system("python " + codedir + "/expression.csv.json.py " + dir)
os.system("python " + codedir + "/cells.csv.py " + dir)
if (expressionRun):
	os.system("python " + codedir + "/expression.csv.py " + dir + " " + size)
