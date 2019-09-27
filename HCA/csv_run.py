import sys, os

def parse(configfile):
	fin = open(configfile, 'r')
	dic ={}
	for line in fin.readlines():
		key, value = line[:-1].split('\t')
		dic[key] = value
	fin.close()
	return dic

if len(sys.argv[:])!= 2:
        print ("python csv_run.py dataDir")
        sys.exit()
    
dir = sys.argv[1]
configfile = dir + '/config'

if not os.path.exists(configfile):
        print ("not config file")
 	sys.exit()

metaDic = parse(configfile)
print (metaDic)
