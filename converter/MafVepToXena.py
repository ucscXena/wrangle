import sys, os, json
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = 'vep.maf$'
columns = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
	'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short', 't_depth', 't_alt_count']

def doChr (data, columnPos): return data[columnPos['Chromosome']]
def doStart (data, columnPos): return data[columnPos['Start_Position']]
def doEnd(data, columnPos): return data[columnPos['End_Position']]
def doRef (data, columnPos): return data[columnPos['Reference_Allele']]
def doAlt (data, columnPos): return data[columnPos['Tumor_Seq_Allele2']]
def doGene (data, columnPos): return data[columnPos['Hugo_Symbol']]
def doEffect (data, columnPos): return data[columnPos['Variant_Classification']]
def doAA (data, columnPos): return data[columnPos['HGVSp_Short']]
def doDNA_VAF (data, columnPos): return str(float(data[columnPos['t_alt_count']]) / float(data[columnPos['t_depth']]))

columnfunctions = {
	'chr': [doChr, 0],
	'start': [doStart, 1],
	'end': [doEnd, 2],
	'reference': [doRef, 3],
	'alt': [doAlt, 4], 
	'gene': [doGene, 5],
	'effect': [doEffect, 6],
	'Amino_Acid_Change': [doAA, 7],
	'DNA_VAF': [doDNA_VAF, 8]
}

def buildjson(assembly, cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='mutationVector'
	J['dataSubtype'] = 'somatic mutation (SNP and INDEL)'
	J['assembly'] = assembly
	J['cohort'] = cohort
	json.dump(J, fout, indent = 4)
	fout.close()

if len(sys.argv[:]) < 5:
	print ("python MafVepToXena.py inputFileDir output cohort \
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

