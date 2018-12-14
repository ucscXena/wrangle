import string, sys
import json
import scipy.stats
import math

# http://onlinestatbook.com/2/chi_square/one-way.html
# chi2: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2.html

def parseObserved(file):
	fin = open(file, 'U')
	dataDic ={}
	cohorts = string.split(fin.readline()[:-1], '\t')[2:] # cohort actually starting on index 2

	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		pathway = data[0]
		if pathway not in dataDic:
			dataDic[pathway] = {}
		else:
			print "Error"
			sys.exit()
		for i in range (2, len(data)):
			value = float(data[i])
			cohort = cohorts[i-2]
			if cohort not in dataDic[pathway]:
				dataDic[pathway][cohort] = value
			else:
				print "Error"
				sys.exit()
	fin.close()
	return dataDic, cohorts

def parseExpected(file):
	fin = open(file, 'U')
	dataExpectedDic ={}
	dataTotalDic = {}
	cohorts = string.split(fin.readline()[:-1], '\t') # cohort actually starting on index 2

	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		pathway = data[0]
		if pathway not in dataExpectedDic:
			dataExpectedDic[pathway] = {}
		else:
			print "Error"
			sys.exit()
		for i in range (2, len(data)):
			E, N = string.split(data[i],',')
			E = float(E)
			N = float(N)
			cohort = cohorts[i]
			if cohort not in dataTotalDic:
				dataTotalDic[cohort] = N
			if cohort not in dataExpectedDic[pathway]:
				dataExpectedDic[pathway][cohort] = E
			else:
				print "Error"
				sys.exit()
	fin.close()
	return dataExpectedDic, dataTotalDic

def chiSqure(expectedDataDic, observedDataDic, totalDataDic, pathways, cohorts, output):
	fout = open(output, 'w')
	fout.write('pathway_label\tpathway_gene_N')
	for cohort in cohorts:
		list = ['cohort', 'Expected', 'Observed', 'N', '%_samples', '%_samples_perGene', 'log2_ObsDivExp',
			'chi_squre_value', '1-prob', 'prob', '-log10_prob']
		list = map(lambda x: cohort + "_" + x, list)
		fout.write('\t'+ string.join(map(str, list), '\t'))
	fout.write('\n')

	for pathway in pathways:
		pathway_label = pathway['golabel']
		pathway_gene_N = len(pathway['gene'])

		fout.write(pathway_label + '\t' + str(pathway_gene_N)) 
		for cohort in cohorts:
			N = totalDataDic[cohort]
			E1 = expectedDataDic[pathway_label][cohort]
			O1 = observedDataDic[pathway_label][cohort]
			if O1 - E1 >=0:
				direction = 1
			else:
				direction = -1
			percent_samples = O1 / N
			percent_samples_pergene = O1 / N / pathway_gene_N
			if O1 >0:
				log2_ObsDivExp = math.log(O1/E1, 2)
			else:
				log2_ObsDivExp = -4 #2E-4
			E2 = N - E1
			O2 = N - O1
			chi_squre_value = (E1 -O1) * (E1 -O1) / E1 + (E2 -O2) * (E2 -O2) / E2
			df = 1
			cdf = scipy.stats.chi2.cdf(chi_squre_value, df)
			prob = 1 - cdf
			if prob != 0.0:
				log10_prob = math.log10(prob)
			else:
				log10_prob = -16 #float('inf')
			list = [cohort, E1, O1, N, percent_samples, percent_samples_pergene, log2_ObsDivExp,
				direction * chi_squre_value, direction * cdf, direction * prob, direction*(-log10_prob)]
			fout.write('\t'+ string.join(map(str, list), '\t'))
		fout.write('\n')
	fout.close()

if len(sys.argv[:]) != 5 and (__name__ == "__main__"):
	print "python chiSqure.py pathway_input expected_input observed_input output"
	print
	sys.exit()

pathway_input  = sys.argv[1]
expectedFile = sys.argv[2]
observedFile = sys.argv[3]
output = sys.argv[4]

# read pathway information
fpathway = open(pathway_input, 'U')
pathways = json.loads(fpathway.read())
fpathway.close()

expectedDataDic, totalDataDic = parseExpected(expectedFile)
observedDataDic, cohorts = parseObserved(observedFile)
chiSqure(expectedDataDic, observedDataDic, totalDataDic, pathways, cohorts, output)
