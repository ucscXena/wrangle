import string, sys

transcriptProbeMap = "gencode.v23.basic.annotation.transcript.probemap"

def processTranscriptProbeMap (transcriptProbeMap):
	fin = open(transcriptProbeMap,'U')
	fin.readline()
	TTSDic = {} # key: hugo values: TTS list
	TTSchrom = {}
	chromGeneList ={} #key: chrom values: gene list
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		gene = data[1]
		strand = data[5]
		chrom = data[2]
		if strand == "+":
			TTS = int(data[3])
		else:
			TTS = int(data[4])
		if gene not in TTSDic:
			TTSDic[gene] =[]
		TTSDic[gene].append(TTS)
		TTSchrom[gene] = chrom

		if chrom not in chromGeneList:
			chromGeneList[chrom] =[]
		if gene not in chromGeneList[chrom]:
			chromGeneList[chrom].append(gene)
	fin.close()
	return [TTSDic, TTSchrom, chromGeneList]

def processPeakInfo(peakinfoFile, TTSDic, TTSchrom, chromGeneList, probeMap):
	fin = open(peakinfoFile,'U')
	fin.readline()
	fout = open(probeMap ,'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	count = 0
	for line in fin.readlines():
		count = count +1
		if count % 20 == 0:
			print count
		data = string.split(line[:-1], '\t')
		chrom = data[0]
		start = int(float(data[1]))
		end = int(float(data[2]))
		id = data[5]
		strand = data[4]
		peak_pos = start + 250
		gene_list =[]
		for gene in chromGeneList[chrom]:
			if TTSchrom[gene] != chrom:
				continue
			for pos in TTSDic[gene]:
				if abs(peak_pos - pos) <= 500 * 1000:
					if gene not in gene_list:
						gene_list.append(gene)
		fout.write(string.join([id, string.join(gene_list,','), chrom, str(start), str(end), strand ],'\t')+'\n')

	fin.close()
	fout.close()


if len(sys.argv[:])!=3:
	print "python buildAllProbeMap.py peak_coord_info_file_in probeMap_file_out"
	print
	sys.exit()
	
peakinfoFile = sys.argv[1]
probeMap = sys.argv[2]

TTSDic, TTSchrom, chromGeneList = processTranscriptProbeMap (transcriptProbeMap)
processPeakInfo(peakinfoFile, TTSDic, TTSchrom, chromGeneList, probeMap)