import string, sys, os
import json

import perSampleEventN

if len(sys.argv[:]) < 5 and (__name__ == "__main__"):
	print "python observed.py pathway_input mode cohort dataInputfile(s)"
	print
	sys.exit()

pathway_input  = sys.argv[1]
mode = sys.argv[2]
cohort = sys.argv[3]
dataInputfiles = sys.argv[4:]

if not os.path.exists(mode):
	os.system("mkdir " + mode)

summaryfile = mode + '/' + cohort +"_pathway_observed"

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