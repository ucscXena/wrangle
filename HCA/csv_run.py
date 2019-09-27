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
    print ("python csv_run.py dataDir")
    sys.exit()

configfile = dir + '/config'

 if not os.exits(configfile):
 	print ("not config file")
 	sys.exit()

metaDic = parse(configfile)
print metaDic