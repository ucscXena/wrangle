import string, sys
import json

import perSampleEventN

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

if len(sys.argv[:]) !=5 and (__name__ == "__main__"):
	print "python  expectedHypergeometric.py pathway_input cohort genome_total_gene_N dataInputfile"
	print
	sys.exit()


pathway_input  = sys.argv[1]
cohort = sys.argv[2]
genome_total_gene_N = int(sys.argv[3])
dataInputfile = sys.argv[4]

perSampleN_input = cohort + "_sampleEvent"
outputfile = cohort +"_sample_expected"
summaryfile = cohort +"_pathway_expected"

#calculate perSample Event number using module
dataDic, genes, genome_total_gene_N, samples = perSampleEventN.perSampleEventN(dataInputfile, genome_total_gene_N)
perSampleEventN.output(perSampleN_input, dataDic, genome_total_gene_N)

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

fperSampleN = open(perSampleN_input, 'U')
fperSampleN.readline()

for line in fperSampleN.readlines():
	sample, N, K = string.split(line[:-1], '\t')
	N = int(N)
	K = int(K)
	fout.write(sample)
	for pathway in pathways:
		prob = 1.0
		n = len(pathway['gene'])
		for i in range (0, n):
			prob = prob * (N-K-i) / (N -i)
		prob = 1 - prob ########### prob (k>=1)
		fout.write('\t'+ str(float(prob)))
		dataDic[pathway['golabel']]= dataDic[pathway['golabel']] + prob
	fout.write('\n')
fperSampleN.close()
fout.close()

fsummaryout = open(summaryfile, 'w')
fsummaryout.write('Expexted\tgene_n\t' + cohort +'\n')
for pathway in pathways:
	fsummaryout.write(pathway['golabel']+'\t' + str(len(pathway['gene'])) + '\t' + str(dataDic[pathway['golabel']]) + ',' + str(len(samples))+'\n')
fsummaryout.close()


