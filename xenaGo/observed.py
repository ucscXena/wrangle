import string, sys
import json

import perSampleEventN

if len(sys.argv[:]) < 4 and (__name__ == "__main__"):
	print "python observed.py pathway_input cohort dataInputfile(s)"
	print
	sys.exit()

pathway_input  = sys.argv[1]
cohort = sys.argv[2]
dataInputfiles = sys.argv[3:]

summaryfile = cohort +"_pathway_observed"

# read pathway information
fpathway = open(pathway_input, 'U')
pathways = json.loads(fpathway.read())
fpathway.close()

fsummaryout = open(summaryfile, 'w')
fsummaryout.write('Observed\tgene_n\t' + cohort +'\n')
for pathway in pathways:
	print pathway['golabel']
	mergedDataDic = perSampleEventN.perSampleEventInPathway(dataInputfiles, pathway['gene'])
	samplesHitN = len(filter(lambda x: x >0, mergedDataDic.values()))
	fsummaryout.write(pathway['golabel']+'\t' + str(len(pathway['gene'])) + '\t' + str(samplesHitN)+'\n')
fsummaryout.close()