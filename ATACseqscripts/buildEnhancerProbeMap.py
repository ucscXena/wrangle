import string, sys

#### gene enhancer linkage for pan-can
enhancerlinkageFile = "ATAC_Peak_Linkage.txt"
probeFile = "TCGA_ATAC_peak.probeMap"
newProbeMap = "TCGA_ATAC_peak.Linkage.probeMap"

#### gene enhancer linkage for brca
#enhancerlinkageFile = "BRCA_ATAC_Peak_Linage.txt"
#probeFile = "brca_peak.probeMap"
#newProbeMap = "brca_enhencer.probeMap"

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
	fin.readline()
	newProbeMapDic = {}
	for line in fin.readlines():
		data = string.split(line[:-1], '\t')
		chrom = data[0]
		start = str(int(float(data[1])))
		end = str(int(float(data[2])))
		targetGene = data[6]
		coord_id = string.join([chrom, start, end],'_')
		try:
			peakID = probeMapDic[coord_id]
			if peakID not in newProbeMapDic:
				newProbeMapDic[peakID] = [[targetGene], chrom, start, end]
			else:
				if targetGene not in newProbeMapDic[peakID][0]:
					newProbeMapDic[peakID][0].append(targetGene)
		except KeyError:
			print coord_id
			continue

	fin.close()

	fout = open(newProbeMap, 'w')
	fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
	for peakID in newProbeMapDic:
		targetGeneS, chrom, start, end = newProbeMapDic[peakID]
		fout.write(string.join([peakID, string.join(targetGeneS, ','), chrom, start, end, '+'],'\t') + '\n')
	fout.close()

probeMapDic = processLowQprobeMap(probeFile)
processLinkagefile(enhancerlinkageFile, probeMapDic)