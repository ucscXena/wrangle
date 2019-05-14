import string, sys, numpy
from sets import Set
import math

def parseZdata(zdata, genes):
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
		gene = line[:string.find(line, '\t')]
		if gene not in genes:
			continue

		data = string.split(line[:-1],'\t')
		zvalues = map(lambda x: float(x), data[1:])
		for i in range (0, N-1):
			if zvalues[i]>3.0:
				zvalues[i] = 3.0
			if zvalues[i] <-3.0:
				zvalues[i] = -3.0
		zDic[gene]= zvalues
	
	fin.close()
	return N-1, samples, zDic

def parseSignatures(signaturefile):
	fin = open(signaturefile,'r')

	sigDic = {}
	name =''
	for line in fin.readlines():
		if line[0]== "#" : #signature name
			name = string.strip(line[1:])
		else:
			genes = string.split(line[1:])
			sigDic[name]= genes
	fin.close()
	return sigDic


def calculateScore(genes,zDic, N):
	total =[]
	for i in range (0, N):
		total.append(0.0)
	count = 0
	for gene in genes:
		if gene not in zDic:
			continue
		count = count + 1
		for i in range (0, N):
			total[i] = total[i] + zDic[gene][i]
	average = map(lambda x : x/count, total)
	return average

def calculateSigScores(sigDic, zDic, N):	
	averageDic = {}
	for signature in sigDic:
		name = signature
		genes = sigDic[name]
		average = calculateScore(genes,zDic, N)
		averageDic[name] = average
	return averageDic

def classification(samples, averageDic):
	# method:  use highest score, and the highest score must > 0.4
	index = averageDic.keys()
	classDic = {}
	classDic['classification']=[]
	for i in range(0, len(samples)):
		sample = samples[i]
		classification = ''
		maxScore = float("-inf")
		for key in index:
			score = averageDic[key][i]
			if score > maxScore:
				maxScore = score
				classification = key
		if maxScore	< 0.4:
			classification = ''
		classDic['classification'].append(classification)
	return classDic

def output_T (samples, dic, outputfile):
	fout = open(outputfile, 'w')
	index = dic.keys()
	fout.write('cell\t' + string.join(index,'\t') +'\n')
	for i in range(0, len(samples)):
		sample = samples[i]
		fout.write(sample)
		for key in index:
			try:
				fout.write('\t' + numpy.format_float_positional(dic[key][i]))
			except:
				fout.write('\t' + dic[key][i])
		fout.write('\n')
	fout.close()

	
if len(sys.argv[:])!= 5:
	print "python signatureScore.py z_input sinaturefile output_score output_classification"
	print
	sys.exit()

zdata = sys.argv[1]
signaturefile = sys.argv[2]
outputfile_score = sys.argv[3]
outputfile_class = sys.argv[4]


sigDic = parseSignatures(signaturefile)

genes =Set([])
for list in sigDic.values():
	genes = genes.union(list)

N, samples, zDic = parseZdata(zdata, genes)
averageDic = calculateSigScores(sigDic, zDic, N)

output_T(samples, averageDic, outputfile_score)

# if this is multi-signature profile, do classification based on signature score 
if len(sigDic.keys()) > 1:
	classDic = classification(samples, averageDic)
	output_T(samples, classDic, outputfile_class)

