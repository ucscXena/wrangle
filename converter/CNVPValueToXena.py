import sys, os, json
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = 'CNVs.p.value.txt$'
columns = ['chr', 'start', 'end', 'copy number']

def doChr (data, columnPos): return data[columnPos['chr']]
def doStart (data, columnPos): return data[columnPos['start']]
def doEnd(data, columnPos): return data[columnPos['end']]
def doCNV (data, columnPos): return data[columnPos['copy number']]

columnfunctions = {
	'chr': [doChr, 0],
	'start': [doStart, 1],
	'end': [doEnd, 2],
	'value': [doCNV, 3]
}

def buildjson(assembly, cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='genomicSegment'
	J['dataSubtype'] = 'copy number segment'
	J['assembly'] = assembly
	J['cohort'] = cohort
	J['colNormalization'] = 'normal2'
	J['max'] = 4
	J['origin'] = 2
	J['thresh'] = 0 
	J["unit"] = "integer copy number",
	json.dump(J, fout, indent = 4)
	fout.close()

if len(sys.argv[:]) < 5:
	print ("python CNVPValueToXena.py inputFileDir output cohort \
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

filearray_convert.fileArrayMafToXena(suffix, columns, inputdir, output, columnfunctions, sample_mapping_file)

buildjson(assembly, cohort, output)