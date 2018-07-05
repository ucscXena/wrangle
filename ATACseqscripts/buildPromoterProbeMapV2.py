import string, sys

def processPeakInfo(peakinfoFile, probeMap):
	fin = open(peakinfoFile,'U')
	fin.readline()
	fout = open(probeMap ,'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")

	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		type = data[7]
		chrom = data[0]
		start = int(float(data[1]))
		end = int(float(data[2]))
		id = data[5]
		strand = data[4]
		
		if type == "Promoter":
			gene = data[11]
		else:
			gene =''
		
		fout.write(string.join([id, gene, chrom, str(start), str(end), strand ],'\t')+'\n')
	fin.close()
	fout.close()


if len(sys.argv[:])!=3:
	print "python buildPromoterProbeMapV2.py peak_coord_info_file_in probeMap_file_out"
	print
	sys.exit()
	
peakinfoFile = sys.argv[1]
probeMap = sys.argv[2]

processPeakInfo(peakinfoFile, probeMap)