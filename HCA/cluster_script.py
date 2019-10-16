import sys,os

def info (small_old_file):
	fin = open(small_old_file,'r')
	mapDic ={}
	for line in fin.readlines():
		data = line[:-1].split('\t')
		xenaSampleID = data[0]
		mapDic [xenaSampleID] = data[1:]
	fin.close()
	return mapDic

def blah (big_mapping_file, smallfile, outputFile, mapDic):
	fin = open(big_mapping_file,'r')
	fout = open(outputFile,'w')

    #header use small file
	fsmall = open(smallfile, 'r')
	line = fsmall.readline()
	N = len(line.split('\t')) - 1
	fout.write(line)
	fsmall.close()

    #data
	while 1:
		line = fin.readline()
		if line =="":
			break
		data = line[:-1].split()
		old_id = data[1]
		new_id = data[0]                
		if old_id not in mapDic:
			fout.write(new_id + '\t' + '\t'* N +'\n')
		else:
			fout.write(new_id + '\t' + '\t'.join(mapDic[old_id]) + '\n')

	fin.close()
	fout.close()

if len(sys.argv[:]) != 4:
	print ("python cluster_script.py cluster_results_small_old_id cluster_results_big_new_id big_mapping(new old)")
	sys.exit()

smallfile = sys.argv[1] # small (old)
output = sys.argv[2] # big (new)
big_mapping_file = sys.argv[3] # big (New old)

mapDic = info(smallfile)
blah (big_mapping_file, smallfile, output, mapDic)


