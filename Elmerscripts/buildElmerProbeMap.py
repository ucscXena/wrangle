import string, sys

illuminProbemapfile ="illuminaMethyl450_hg38_GDC"

def probeMapParse(illuminProbemapfile):
	probeDic ={}
	fin = open(illuminProbemapfile, 'U')
	fin.readline()
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line[:-1],'\t')
		id = data[0]
		chrom = data[2]
		start = data[3]
		end = data[4]
		probeDic[id] = [chrom, start, end]
	fin.close()
	return probeDic

def elmerFileParse (elemerfile):
	elemerGeneDic ={}
	fin = open(elemerfile, 'U')
	
	for line in fin.readlines():
		data= string.split(line[:-1],',')
		id = data[0]
		if id[:2] != "cg":
			continue
		gene =data[2]
		if id not in elemerGeneDic:
			elemerGeneDic[id]=[]
		if gene not in elemerGeneDic[id]:
			elemerGeneDic[id].append(gene)
	fin.close()
	return elemerGeneDic

def buildElmerProbeMap(probeDic, elemerGeneDic,outputProbeMap):
	fout = open(outputProbeMap, 'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	for id in elemerGeneDic:
		genes = string.join(elemerGeneDic[id],',')
		chrom, start, end = probeDic[id]
		fout.write(string.join([id, genes, chrom, start, end, '.'],'\t')+'\n')
	fout.close()

if len(sys.argv[:]) != 3:
	print "python buildElmerProbeMap.py elemer_filein probeMap_fileout"
	print
	sys.exit()

elemerfile = sys.argv[1]
outputProbeMap= sys.argv[2]

probeDic = probeMapParse(illuminProbemapfile)
elemerGeneDic = elmerFileParse (elemerfile)

buildElmerProbeMap(probeDic, elemerGeneDic,outputProbeMap)