import string, copy, numpy
import re

#zdata = 'z.txt'
zdata = 'infercnv.observations.txt'
probeMap = 'NI03_CNV_hg19_genes_ordered_correct_noXY.probeMap'
cnvOutput = 'cnv.txt'

window =100  			# gene number window for averaging
jump = 10 				# select every N genes to output , display limitation
centerCell = False		# center cell value

def parseZdata(zdata):
	fin =open(zdata,'r')

	line = fin.readline()
	data = string.split(line[:-1],'\t')
	N = len(data)
	samples = data[1:]

	zDic ={}
	
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line[:-1],'\t')
		gene = data[0]
		zvalues = map(lambda x: float(x), data[1:])
		for i in range (0, N-1):
			if zvalues[i]>3.0:
				zvalues[i] = 3.0
			if zvalues[i] <-3.0:
				zvalues[i] = -3.0
		zDic[gene]= zvalues
	
	fin.close()
	return N-1, samples, zDic


def parseProbeMap(probeMap, zDic):
	fin =open(probeMap,'r')
	fin.readline()
	genePos = {}
	for line in fin.readlines():
		data = string.split(line[:-1],'\t')
		gene = data[1]
		if gene not in zDic:
			continue
		chr = data[2]
		if string.find(chr,'_') != -1:
			continue
		start = float(data[3])

		# convert chr1 CHR1 to 1, numbers only
		chr = re.sub(re.escape('chr'), '', chr, flags=re.IGNORECASE)
		try:
			chr = int(chr) # ignore X Y
		except:
			pass
		if chr not in genePos:
			genePos[chr] =[]
		genePos[chr].append([start, gene])
	fin.close()

	for chr in genePos:
		genePos[chr].sort()
	return genePos

def calculateCNVtotal(genes, zDic, N):
	total =[]
	for i in range (0, N):
		total.append(0.0)
	for gene in genes:
		if gene not in zDic:
			continue
		for i in range (0, N):
			total[i] = total[i] + zDic[gene][i]
	return total

def calculateCNVincrement(genesAdd, genesMinus, zDic, N, total):
	for i in range (0, N):
		for j in range(0, len(genesAdd)):
			total[i] = total[i] + zDic[genesAdd[j]][i] - zDic[genesMinus[j]][i]
	return total

def convertExpToCNV(zDic, genePos, N):
	cnvDic = {}
	cnvDicKeyOrder = [] # because we compute from chr1,chr2, we want to keep the order
	total =[]
	for i in range (0, N):
		total.append(0.0)

	chromOrder = genePos.keys.sort()
	for chr in chromOrder:
		genes = map(lambda x: x[1], genePos[chr][:window])
		total = calculateCNVtotal(genes, zDic, N)
		
		i = 0
		cnvDic[chr + '_' + str(i)] = total
		
		while 1:
			i = i + jump
			if i >  len(genePos[chr]) - window:
				break
			genesAdd = map(lambda x: x[1], genePos[chr][i + window -1 - jump : i + window - 1])
			genesMinus = map(lambda x: x[1], genePos[chr][i-jump: i])
			
			total = calculateCNVincrement(genesAdd, genesMinus, zDic, N, total)
			probe = 'chr' + str(chr) + '_' + str(i)
			cnvDic[probe] = copy.deepcopy(total)
			cnvDicKeyOrder.append(probe)

#		for i in range (jump, len(genePos[chr]) - window, jump):
#			genesAdd = map(lambda x: x[1], genePos[chr][i + window -1 : i + window - 1 + jump])
#			genesMinus = map(lambda x: x[1], genePos[chr][i-jump: i])

#			total = calculateCNVincrement(genesAdd, genesMinus, zDic, N, total)
#			cnvDic[chr + '_' + str(i)] = copy.deepcopy(total)
	return cnvDic, cnvDicKeyOrder

def centerPerCell(cnvDic, N):
	posList = cnvDic.keys()
	for i in range(0, N):
		total = 0.0
		for pos in posList:
			total = total + cnvDic[pos][i]
		average = total / len(posList)
		for pos in posList:
			cnvDic[pos][i] = cnvDic[pos][i] - average
	return cnvDic

def output (samples, cnvDic, cnvOutput, cnvDicKeyOrder):
	fout = open(cnvOutput, 'w')
	fout.write('cell\t' + string.join(samples,'\t') +'\n')
	for key in cnvDicKeyOrder:
		fout.write(key+'\t')
		fout.write(string.join(map(lambda x: numpy.format_float_positional(x/window), cnvDic[key]),'\t'))
		fout.write('\n')
	fout.close()

N, samples, zDic = parseZdata(zdata)
genePos = parseProbeMap(probeMap, zDic)
cnvDic, cnvDicKeyOrder = convertExpToCNV(zDic, genePos, N)
if centerCell:
	cnvDic = centerPerCell(cnvDic, N)
output (samples, cnvDic, cnvOutput, cnvDicKeyOrder)
