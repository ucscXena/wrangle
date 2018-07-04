import string, sys

#### gene enhancer linkage for pan-can
# enhancerlinkageFile = "ATAC_Peak_Linkage.txt"
# probeFile = "TCGA_ATAC_peak.probeMap"
# newProbeMap = "TCGA_ATAC_peak.Linkage.probeMap"

#### gene enhancer linkage for brca
enhancerlinkageFile = "BRCA_ATAC_Peak_Linage.txt"
probeFile = "brca_peak.probeMap"
newProbeMap = "brca_enhencer.probeMap"

def processLowQprobeMap(probeFile):
	fin = open(probeFile,'U')
	fin.readline()
	probeMapDic ={} #key: chr_start_end value: probeID
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		peakID = data[0]
		chrom = data[2]
		start = data[3]
		end = data[4]
		probeMapDic[string.join([chrom, start, end],'_')] = peakID
	fin.close()
	return probeMapDic

def processLinkagefile(enhancerlinkageFile, probeMapDic):
	fin = open(enhancerlinkageFile,'U')
	fout = open(newProbeMap, 'w')
	fin.readline()
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		enhancerID = data[-1]
		chrom = data[0]
		start = str(int(float(data[1])))
		end = str(int(float(data[2])))
		targetGene = data[6]
		coord_id = string.join([chrom, start, end],'_')
		try:
			peakID = probeMapDic[coord_id]
			fout.write(string.join([peakID, targetGene, chrom, start, end, '+'],'\t') + '\n')
		except KeyError:
			print coord_id
	fin.close()
	fout.close()

probeMapDic = processLowQprobeMap(probeFile)
processLinkagefile(enhancerlinkageFile, probeMapDic)