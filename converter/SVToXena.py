import sys, os, json
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = '.vcf$'
columns = ['CHROM', 'POS', 'REF', 'ALT', 'INFO']

def parseInfo(line):
	dic ={}
	for item in line.split(";"):
		data = item.split("=")
		if len(data) == 2:
			key, value = 
		dic[key] = value
	return dic

def doChr (data, columnPos): return data[columnPos['CHROM']]
def doStart (data, columnPos): return data[columnPos['POS']]
def doAlt (data, columnPos): return data[columnPos['ALT']]
def doEffect (data, columnPos): return parseInfo(data[columnPos['INFO']])['SVTYPE']
def doEnd(data, columnPos): 
	if parseInfo(data[columnPos['INFO']])['SVTYPE'] == 'BND':
		return data[columnPos['POS']]
	else:
		return parseInfo(data[columnPos['INFO']])['END']

def doRef (data, columnPos):
	if parseInfo(data[columnPos['INFO']])['SVTYPE'] == 'BND':
		return data[columnPos['REF']]
	else:
		return data[columnPos['POS']] + '-' + parseInfo(data[columnPos['INFO']])['END']



columnfunctions = {
	'chr': [doChr, 0],
	'start': [doStart, 1],
	'end': [doEnd, 2],
	'reference': [doRef, 3],
	'alt': [doAlt, 4],
	'effect': [doEffect, 5]
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