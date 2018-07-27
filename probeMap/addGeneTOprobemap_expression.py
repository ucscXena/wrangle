import string, sys

#overlapping exon region only, strand specific

def processTranscriptProbeMap (exonProbeMap):
	fin = open(transcriptProbeMap,'U')
	fin.readline()
	regionDic = {} # key: hugo values: region list
	chromDic = {}
	strandDic= {}
	chromGeneList ={} #key: chrom values: gene list
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		gene = data[1]
		strand = data[5]
		chrom = data[2]
		start = int(data[3])
		end = int(data[4])
		strandDic[gene] = strand
		chromDic[gene] = chrom
		if gene not in regionDic:
			regionDic[gene] =[]
		regionDic[gene].append([start,end])

		if chrom not in chromGeneList:
			chromGeneList[chrom] =[]
		if gene not in chromGeneList:
			chromGeneList[chrom].append(gene)

	fin.close()
	return [regionDic, strandDic, chromDic, chromGeneList]

def augmentProbeMap(probeMap_coord, regionDic, strandDic, chromDic, chromGeneList, probeMap_out):
	fin = open(probeMap_coord,'U')
	fin.readline()
	fout = open(probeMap_out ,'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	count = 0
	for line in fin.readlines():
		count = count +1
		if count % 20 == 0:
			print count
		data = string.split(line[:-1], '\t')
		chrom = data[2]
		start = int(data[3])
		end = int(data[4])
		id = data[0]
		strand = data[5]

		gene_list =[]
		for gene in chromGeneList[chrom]:
			if chromDic[gene] != chrom:
				continue
			if strandDic[gene] != strand:
				continue
			for exon in regionDic[gene]:
				eStart, eEnd = exon
				if start < eEnd and end > eStart:
					gene_list.append(gene)
					break
		fout.write(string.join([id, string.join(gene_list,','), chrom, str(start), str(end), strand ],'\t')+'\n')

	fin.close()
	fout.close()


if len(sys.argv[:])!=4:
	print "python addGeneToprobemap_expression.py exon_annotation(probeMap) probeMap_in probeMap_out"
	print "#overlapping exon region only, strand specific"
	print
	sys.exit()

exonProbeMap = sys.argv[1]
probeMap_coord = sys.argv[2]
probeMap_out = sys.argv[3]

regionDic, strandDic, chromDic, chromGeneList = processTranscriptProbeMap (exonProbeMap)
augmentProbeMap(probeMap_coord, regionDic, strandDic, chromDic, chromGeneList, probeMap_out)