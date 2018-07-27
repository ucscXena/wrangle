import string, sys

if len(sys.argv[:]) !=3:
	print "python xenaGeneModelToProbeMap.py xeneGeneGood_filein hugo_probemap_fileout"
	print
	sys.exit()

goodGeneFile = open(sys.argv[1],'U')
goodGeneFile.readline()

probeMapfile = open(sys.argv[2],'w')
probeMapfile.write(string.join(['id', 'gene', 'chrom', 'chromStart', 'chromEnd', 'strand', 'thickStart', 'thickEnd', 'blockCount', 'blockSizes', 'blockStarts'], '\t') + '\n')


for line in goodGeneFile.readlines():
	data = string.split(line[:-1],'\t')
	gene = data[12]
	chrom = data[2]
	strand = data[3]
	start = int(data[4])
	chromStart = str(start + 1)
	chromEnd = data[5]
	thickStart = str(int(data[6]) +1)
	thickEnd = data[7]
	blockCount = data[8]

	exonStarts = map(lambda x : int(x), string.split(data[9],',')[:-1]) 
	exonEnds = map(lambda x : int(x), string.split(data[10],',')[:-1])

	assert int(blockCount) == len(exonStarts)	
	assert int(blockCount) == len(exonEnds)

	blockSizes = map(lambda x: x[1] - x[0] , zip(exonStarts, exonEnds))
	blockSizes = string.join(map(lambda x: str(x) ,blockSizes),',')

	blockStarts = map(lambda x: x - start, exonStarts)
	blockStarts = string.join(map(lambda x: str(x) ,blockStarts),',')

	probeMapfile.write(string.join([gene, gene, chrom , chromStart, chromEnd, strand, thickStart, thickEnd, blockCount, blockSizes, blockStarts], '\t') + '\n')

goodGeneFile.close()
probeMapfile.close()
