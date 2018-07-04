import string, sys

transcriptProbeMap = "gencode.v23.basic.annotation.transcript.probeMap"

def processTranscriptProbeMap (transcriptProbeMap):
	fin = open(transcriptProbeMap,'U')
	fin.readline()
	TTSDic = {} # key: hugo values: TTS list
	TTSchrom = {}
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		gene = data[2]
		strand = data[5]
		chrom = data[1]
		if strand == "+":
			TTS = int(data[3])
		else:
			TTS = int(data[4])
		if gene not in TTSDic:
			TTSDic[gene] =[]
		TTSDic[gene].append(TTS)
		TTSchrom[gene] = chrom
	fin.close()
	return [TTSDic, TTSchrom]

def processPeakInfo(peakinfoFile, TTSDic, TTSchrom, probeMap):
	fin = open(peakinfoFile,'U')
	fin.readline()
	fout = open(probeMap ,'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	count = 0
	for line in fin.readlines():
		count = count +1
		print count
		data = string.split(line[:-1], '\t')
		chrom = data[0]
		start = int(float(data[1]))
		end = int(float(data[2]))
		id = data[5]
		strand = data[4]
		peak_pos = start + 250
		gene_list =[]
		for gene in TTSchrom:
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

TTSDic, TTSchrom = processTranscriptProbeMap (transcriptProbeMap)
processPeakInfo(peakinfoFile, TTSDic, TTSchrom, probeMap)