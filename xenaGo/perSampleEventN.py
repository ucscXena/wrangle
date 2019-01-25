import string, sys
from sets import Set

def PositiveGisticThreshold (value):
	if value >=2.0 or value <=-2.0:
		return 1
	return 0

def PositiveProteinMutation (value):
	if value not in ['synonymous_variant',
		'Synonymous Variant',
		'Synonymous',
		'Silent',
		'stop_retained_variant',
		'lincRNA',
		'RNA',
		'exon_variant',
		'upstream_gene_variant',
		'downstream_gene_variant',
		"5'Flank",
		"3'Flank",
		"3'UTR",
		"5'UTR",
		'5_prime_UTR_variant',
		'3_prime_UTR_variant',
		'Complex Substitution',
		'intron_variant',
		'intron',
		'Intron',
		'intergenic_region',
		'IGR']:
		return 1
	return 0

def parse_perSampleEvent(file):
	fin = open(file, 'r')
	header=	fin.readline()	
	
	type = 'GisticThreshold'
	dataDic = {} # key : sample value: total event for that sample
	genes =[]

	data = 	string.split(header[:-1], '\t')
	for value in data[1:]:  # very hacky parser
		if value == 'alt':
			type = 'Mutation'
			break

	# parse GisticThreshold file
	if type == 'GisticThreshold':
		samples = data[1:]
		N = len(samples) # total number of samples
		totalEvent = []

		for i in range (0, N):
			totalEvent.append(0)

		while 1:
			line = fin.readline()
			if line =='':
				fin.close()
				break
			data = string.split(line[:-1], '\t')
			gene = data[0]
			if gene not in genes:
				genes.append(gene)

			data = map(PositiveGisticThreshold, map(float, data[1:]))
			for i in range(0, N):
				totalEvent[i] = totalEvent[i]+ data[i]
		
		for i in range(0, N):
			sample = samples [i]
			totalE = totalEvent[i]
			dataDic[sample] = totalE

	# parse xena mutation vector file
	elif type == 'Mutation':
		headerDic ={}
		headerDic['sample'] = 0	 	
		samples = []
		for i in range (1, len(data)):
			if data[i] == 'effect':
				headerDic['effect'] = i
			if data[i] == 'gene':
				headerDic['gene'] = i
		if 'effect' not in headerDic or 'gene' not in headerDic:
			print 'ERROR: missing header: effect or gene'
			sys.exit()
		while 1:
			line = fin.readline()
			if line =='':
				fin.close()
				break
			data = string.split(line[:-1], '\t')
			sample = data[0]
			if sample not in samples:
				samples.append(sample)
				dataDic[sample]={}
			effect = data[headerDic['effect']]
			gene = data[headerDic['gene']]
			if gene not in genes:
				genes.append(gene)
			if PositiveProteinMutation(effect):
				if gene not in dataDic[sample]:
					dataDic[sample][gene] = 1
		for sample in dataDic:
			dataDic[sample] = reduce(lambda x, y : x+y , dataDic[sample].values(), 0)
	
	return type, dataDic, genes, samples


def parse_perSampleEventInPathway(file, genes):
	fin = open(file, 'r')
	header=	fin.readline()	
	
	type = 'GisticThreshold'
	dataDic = {} # key : sample value: total event for that sample in pathway

	data = 	string.split(header[:-1], '\t')
	for value in data[1:]:  # very hacky parser
		if value == 'alt':
			type = 'Mutation'
			break

	# parse GisticThreshold file
	if type == 'GisticThreshold':
		samples = data[1:]
		N = len(samples) # total number of samples
		totalEvent = []

		for i in range (0, N):
			totalEvent.append(0)

		while 1:
			line = fin.readline()
			if line =='':
				fin.close()
				break
			data = string.split(line[:-1], '\t')
			gene = data[0]
			if gene not in genes:
				continue
			data = map(PositiveGisticThreshold, map(float, data[1:]))
			for i in range(0, N):
				totalEvent[i] = totalEvent[i]+ data[i]
		
		for i in range(0, N):
			sample = samples [i]
			totalE = totalEvent[i]
			dataDic[sample] = totalE

	# parse xena mutation vector file
	elif type == 'Mutation':
		headerDic ={}
		headerDic['sample'] = 0	
		samples =[]	
		for i in range (1, len(data)):
			if data[i] == 'effect':
				headerDic['effect'] = i
			if data[i] == 'gene':
				headerDic['gene'] = i
		if 'effect' not in headerDic or 'gene' not in headerDic:
			print 'ERROR: missing header: effect or gene'
			sys.exit()
		while 1:
			line = fin.readline()
			if line =='':
				fin.close()
				break
			data = string.split(line[:-1], '\t')
			sample = data[0]
			if sample not in samples:
				samples.append(sample)
			effect = data[headerDic['effect']]
			gene = data[headerDic['gene']]
			if gene not in genes:
				continue
			if PositiveProteinMutation(effect):
				if sample not in dataDic:
					dataDic[sample]={}
				if gene not in dataDic[sample]:
					dataDic[sample][gene] = 1
		for sample in dataDic:
			dataDic[sample] = reduce(lambda x, y : x+y , dataDic[sample].values())
	
	return type, dataDic, genes, samples


def output(outputfile, dataDic, genome_total_gene_N):
	fout = open(outputfile ,'w')
	fout.write('sample\ttotal_pop_N\tevent_K\n')
	for sample in dataDic:
		fout.write(sample + '\t' + str(genome_total_gene_N) + '\t' + str(dataDic[sample]) + '\n')
	fout.close()

def perSampleEventN (inputfile, genome_total_gene_N):
	type, dataDic, genes, samples = parse_perSampleEvent(inputfile) # key : sample value: total event for that sample
	print inputfile, len(genes), len(samples)

	if genome_total_gene_N == 0:
		genome_total_gene_N = len(genes)

	return dataDic, genes, genome_total_gene_N, samples

def perSampleEventInPathway (inputfiles, pathway_genes):
	mergedDataDic = {}
	mergedSamples = Set ([])
	allDataCollection =[] # list of type, dataDic, genes

	for file in inputfiles:
		type, dataDic, genes, samples = parse_perSampleEventInPathway(file, pathway_genes) # key : sample value: total event for that sample
		allDataCollection.append([type, dataDic, genes])
		if len(mergedSamples)== 0:
			mergedSamples = Set(samples)
		else:
			mergedSamples = mergedSamples.intersection(Set(samples))

	for sample in mergedSamples:
		mergedDataDic[sample] = 0
		for list in allDataCollection:
			type, dataDic, genes = list	
			if sample in dataDic:
				mergedDataDic[sample] = mergedDataDic[sample] + dataDic[sample]
	return mergedDataDic


if len(sys.argv[:]) != 4 and (__name__ == "__main__"):
	print "python  perSampleEventN.py output genome_total_gene_N inputfile"
	print
	sys.exit()

if __name__ == "__main__":
	inputfiles = sys.argv[3]
	genome_total_gene_N = int(sys.argv[2])
	outputfile = sys.argv[1]
	mergedDataDic, mergedGenes, mergedSampleTotalGenes, mergedSamples = perSampleEventN (inputfile, genome_total_gene_N)
	output(outputfile, mergedDataDic, mergedSampleTotalGenes) 
	print "genes:", len(mergedGenes)
