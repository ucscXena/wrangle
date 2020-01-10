import sys, os, json, datetime
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = 'abundance.tsv$'
columns = ["est_counts", "tpm"]
pseudocounts = [1, 0.001]
maxNum = 500

def buildExpJson(cohort, output_prefix, probeMapfilename):
	for i in range (0, len(columns)):
		column = columns[i]
		pseudocount = pseudocounts[i]
		output = output_prefix + '_' + column + '.txt.json'
		fout = open(output, 'w')
		J = {}
		J['type'] ='genomicMatrix'
		if probeMapfilename:
			J[':probeMap'] = probeMapfilename
		J['unit'] = 'log2('+ column + '+' + str(pseudocount) +')'
		J['pseudocount'] = pseudocount
		J['cohort'] = cohort
		J['version'] = datetime.date.today().isoformat()
		J['dataSubtype'] = 'transcript expression'
		json.dump(J, fout, indent = 4)
		fout.close()

if len(sys.argv[:]) < 4:
	print ("python KallistoToXena.py inputFileDir output_prefix cohort \
		probeMapfilename(optional) sample_id_mapping(optional, id file_id)")
	print ()
	sys.exit(1)

inputdir = sys.argv[1]
output_prefix = sys.argv[2]
cohort = sys.argv[3]

if len(sys.argv[:]) > 4:
	probeMapfilename = sys.argv[4]
else:
	probeMapfilename = None

if len(sys.argv[:]) > 5 :
	sample_mapping_file = sys.argv[5]
else:
	sample_mapping_file = None

filearray_convert.fileArrayExpToXena(suffix, columns, maxNum, inputdir, output_prefix, pseudocounts, 
	sample_mapping_file)

buildExpJson(cohort, output_prefix, probeMapfilename)