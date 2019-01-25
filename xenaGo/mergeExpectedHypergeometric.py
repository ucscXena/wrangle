import string, sys
import json

# Hypergeometric distribution https://en.wikipedia.org/wiki/Hypergeometric_distribution
# https://en.wikipedia.org/wiki/Combination 
# prob(k>1)  =  1 - probability of zero event
# N is the population size: 2* gene count
# K is the number of success states in the population: sample event count
# n is the number of draws (i.e. quantity drawn in each trial) : gene number in a pathway
# k is the number of observed successes: 0 
# prob(k=0) = (N-K) * (N-K-1) * ... * (N-K-(n-1))/
#          N * (N-1) * ... * (N- (n-1))
# goal: calculate prob(k>=1)

from itertools import combinations 

def addIndepProb (prob_list):	# p = PA + PB - PAB, etc
	total_prob = 0.0
	for i in range (1, len(prob_list)+1):
		if (i % 2):
			sign = 1.0
		else:
			sign = -1.0
		for combo in combinations(prob_list, i):
			total_prob = total_prob + reduce(lambda x,y: x*y, combo, 1) * sign
	return total_prob


if len(sys.argv[:]) < 4 and (__name__ == "__main__"):
	print "python  mergeExpectedHypergeometric.py pathway_input cohort genome_total_gene_N perSampleNFiles"
	print
	sys.exit()


pathway_input  = sys.argv[1]
cohort = sys.argv[2]
perSampleNInputfiles = sys.argv[3:]

perSampleN_output = cohort + "_sampleEvent"
outputfile = cohort +"_sample_expected"
summaryfile = cohort +"_pathway_expected"

#calculate perSample Event number by merging perSample Event number from multiple filters
mergedDataDic = {}
mergedSampleTotalGenes = {}
mergedSamples ={}

for file in perSampleNInputfiles:
	fperSampleN = open(file, 'U')
	fperSampleN.readline()
	for line in fperSampleN.readlines():
		sample, N, K = string.split(line[:-1], '\t')
		N = int(N)
		K = int(K)

		if sample not in mergedDataDic:
			mergedDataDic[sample] = []
		if sample not in mergedSampleTotalGenes:
			mergedSampleTotalGenes[sample] =[]
		if sample not in mergedSamples:
			mergedSamples[sample] = 0

		mergedDataDic[sample].append(K)
		mergedSampleTotalGenes[sample].append(N)	
		mergedSamples[sample] = mergedSamples[sample] + 1

	fperSampleN.close()

# output merged perSampleN file
fout = open(perSampleN_output, 'w')
fout.write("sample")
for i in range(0,len(perSampleNInputfiles)):
	fout.write("\ttotal_pop_N\tevent_K")
fout.write("\n")

for sample in mergedSamples:
	if mergedSamples[sample] == len(perSampleNInputfiles):   # only intersection samples
		fout.write(sample)
		for i in range(0,len(perSampleNInputfiles)):
			K = mergedDataDic[sample][i]
			N = mergedSampleTotalGenes[sample][i]
			fout.write('\t' + str(N) + '\t' + str(K))
		fout.write("\n")
fout.close()

# read pathway information
fpathway = open(pathway_input, 'U')
pathways = json.loads(fpathway.read())
fpathway.close()

dataDic = {} # key: pathway, value : sum of all samples probability

fout = open(outputfile, 'w')
#fout.write('sample\ttotal_pop_N\tevent_K')
fout.write('sample')
for pathway in pathways:
	pathway_label = pathway['golabel'] + "_" + str(len(pathway['gene']))
	fout.write('\t'+ pathway_label)
	if pathway['golabel'] not in dataDic:
		dataDic[pathway['golabel']] = 0.0
fout.write('\n')

intersection_sample_N = 0
for sample in mergedSamples:
	if mergedSamples[sample] == len(perSampleNInputfiles):   # only intersection samples
		intersection_sample_N = intersection_sample_N + 1

		fout.write(sample)			
		for pathway in pathways:
			sample_prob = []
			for i in range(0,len(perSampleNInputfiles)):
				K = mergedDataDic[sample][i]
				N = mergedSampleTotalGenes[sample][i]	
				prob = 1.0
				n = len(pathway['gene'])
				for i in range (0, n):
					prob = prob * (N-K-i) / (N -i)
				prob = 1 - prob ########### prob (k>=1)
				sample_prob.append(prob)  # combined probability from multiple filters
			total_prob = addIndepProb(sample_prob)
			fout.write('\t'+ str(float(total_prob)))
			dataDic[pathway['golabel']]= dataDic[pathway['golabel']] + total_prob
		fout.write('\n')
fout.close()

fsummaryout = open(summaryfile, 'w')
fsummaryout.write('Expexted\tgene_n\t' + cohort +'\n')
for pathway in pathways:
	fsummaryout.write(pathway['golabel']+'\t' + str(len(pathway['gene'])) + '\t' + str(dataDic[pathway['golabel']]) + ',' + str(intersection_sample_N)+'\n')
fsummaryout.close()

