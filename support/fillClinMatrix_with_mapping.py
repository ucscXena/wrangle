import sys,os

def mapping (mfile):
	fin = open(mfile,'r')
	mapDic ={}
	for line in fin.readlines():
		data = line[:-1].split()
		xenaSampleID = data[0]
		mapDic [xenaSampleID] = line
	fin.close()
	print (mapDic)
	return mapDic

def blah (bigfile, smallfile, outputFile, mapDic):
	fin = open(bigfile,'r')
	fout = open(outputFile,'w')

    #header use small file
	fsmall = open(smallfile, 'r')
	line = fin.readline()
	N = len(line.split('\t')) - 1
	fout.write(line)
	fsmall.close()

    #data
	while 1:
		line = fin.readline()
		if line =="":
			break
		data = line[:-1].split()
		id = data[0]
		if id not in mapDic:
			fout.write(id + '\t' + '\t'* N +'\n')
		else:
			fout.write(mapDic[id])

	fin.close()
	fout.close()

if len(sys.argv[:]) != 4:
	print ("python fillClinMatrix_with_big_useNULL.py input_small(clin, Mutaiton, etc) outupt big(first col id, no header)")
	sys.exit()

smallfile = sys.argv[1] # small
output = sys.argv[2]
bigfile = sys.argv[3] # big

mapDic = mapping(smallfile)
blah (bigfile, smallfile, output, mapDic)


