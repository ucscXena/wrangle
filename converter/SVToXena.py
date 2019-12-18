import sys, os, json
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = '.vcf$'
columns = ['CHROM', 'POS', 'REF', 'ALT']

def doChr (data, columnPos): return data[columnPos['CHROM']]
def doStart (data, columnPos): return data[columnPos['POS']]
def doEnd(data, columnPos): return data[columnPos['POS']]
def doRef (data, columnPos): return data[columnPos['REF']]
def doAlt (data, columnPos): return data[columnPos['ALT']]

columnfunctions = {
	'chr': [doChr, 0],
	'start': [doStart, 1],
	'end': [doEnd, 2],
	'reference': [doRef, 3],
	'alt': [doAlt, 4]
}

def buildjson(assembly, cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='mutationVector'
	J['dataSubtype'] = 'somatic structural variant'
	J['assembly'] = assembly
	J['cohort'] = cohort
	json.dump(J, fout, indent = 4)
	fout.close()

if len(sys.argv[:]) < 5:
	print ("python SVToXena.py inputFileDir output cohort \
		assembly sample_id_mapping(optional, id file_id)")
	print ()
	sys.exit(1)

inputdir = sys.argv[1]
output = sys.argv[2]
cohort = sys.argv[3]
assembly = sys.argv[4]

if len(sys.argv[:]) > 5 :
	sample_mapping_file = sys.argv[5]
else:
	sample_mapping_file = None

filearray_convert.fileArrayVcfToXena(suffix, columns, inputdir, output, columnfunctions, sample_mapping_file)

buildjson(assembly, cohort, output)