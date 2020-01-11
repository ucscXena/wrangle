import sys, os, json, datetime
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = 'arriba.fusions.tsv$'
columns = ['breakpoint1', 'breakpoint2', 
	'gene1', 'gene2', 
	'site1', 'site2', 
	'strand1(gene/fusion)', 'strand2(gene/fusion)',
	'direction1', 'direction2', 
	'type', 'fusion_transcript', 'reading_frame', 'confidence']

compliment = {
	'A' : 'T',
	'T' : 'A',
	'G' : 'C',
	'C' : 'G'
}

def parsePos(genomicPos):
	Chr, pos = genomicPos.split(':')
	if Chr[0:3].upper() != 'CHR':
		Chr = 'chr' + Chr
	return {
		'Chr': Chr,
		'pos': pos
	}

def doChr1 (data, columnPos): return parsePos(data[columnPos['breakpoint1']])['Chr']
def doPos1 (data, columnPos): return parsePos(data[columnPos['breakpoint1']])['pos']
def doRef1 (data, columnPos):
	DNA = data[columnPos['fusion_transcript']].split('|')[0][-1]
	strand = data[columnPos['strand1(gene/fusion)']].split('/')[1]
	if strand == '-':
		DNA = compliment[DNA]
	return DNA

def doAlt1 (data, columnPos):
	ref = doRef1(data, columnPos)
	bp2 = data[columnPos['breakpoint2']]
	d1 = data[columnPos['direction1']]
	d2 = data[columnPos['direction2']]
	if d1 == 'downstream':
		if d2 == 'upstream': # d2 is the downstream piece
			alt = ref + '[' + bp2 + '['
		if d2 == 'downstream': # d2 is the upstream piece
			alt = ref + ']' + bp2 + ']'
	if d1 == 'upstream':
		if d2 == 'downstream' : # d2 is the upstream piece
			alt = ']' + bp2 + ']' + ref
		if d2 == 'upstream': # d2 is the downstream piece
			alt = '[' + bp2 + '[' + ref
	return alt

def doGene1(data, columnPos):
	if data[columnPos['site1']] != 'intergenic':
		return data[columnPos['gene1']]
	else:
		return 'intergenic'

def doAltGene1(data, columnPos):
	if data[columnPos['site2']] != 'intergenic':
		return data[columnPos['gene2']]
	else:
		return 'intergenic'

def doEffect (data, columnPos): return data[columnPos['type']]
def doReadingFrame (data, columnPos): return data[columnPos['reading_frame']]
def doConfidence(data, columnPos): return data[columnPos['confidence']]

def doChr2 (data, columnPos): return parsePos(data[columnPos['breakpoint2']])['Chr']
def doPos2 (data, columnPos): return parsePos(data[columnPos['breakpoint2']])['pos']
def doRef2 (data, columnPos):
	DNA = data[columnPos['fusion_transcript']].split('|')[1][0]
	strand = data[columnPos['strand2(gene/fusion)']].split('/')[1]
	if strand == '-':
		DNA = compliment[DNA]
	return DNA
def doAlt2 (data, columnPos):
	ref = doRef2(data, columnPos)
	bp2 = data[columnPos['breakpoint1']]
	d1 = data[columnPos['direction2']]
	d2 = data[columnPos['direction1']]
	if d1 == 'downstream':
		if d2 == 'upstream': # d2 is the downstream piece
			alt = ref + '[' + bp2 + '['
		if d2 == 'downstream': # d2 is the upstream piece
			alt = ref + ']' + bp2 + ']'
	if d1 == 'upstream':
		if d2 == 'downstream' : # d2 is the upstream piece
			alt = ']' + bp2 + ']' + ref
		if d2 == 'upstream': # d2 is the downstream piece
			alt = '[' + bp2 + '[' + ref
	return alt

def doGene2(data, columnPos):
	if data[columnPos['site2']] != 'intergenic':
		return data[columnPos['gene2']]
	else:
		return 'intergenic'

def doAltGene2(data, columnPos):
	if data[columnPos['site1']] != 'intergenic':
		return data[columnPos['gene1']]
	else:
		return 'intergenic'

columnfunctions = {
	'chr': [doChr1, 0],
	'start': [doPos1, 1],
	'end': [doPos1, 2],
	'reference': [doRef1, 3],
	'alt': [doAlt1, 4],
	'gene': [doGene1, 5],
	'altGene': [doAltGene1, 6],
	'effect': [doEffect, 7],
	'reading_frame': [doReadingFrame, 8],
	'confidence' :[doConfidence, 9]
}

paired_columnfunctions = {
	'chr': [doChr2, 0],
	'start': [doPos2, 1],
	'end': [doPos2, 2],
	'reference': [doRef2, 3],
	'alt': [doAlt2, 4],
	'gene2': [doGene1, 5],
	'altGene2': [doAltGene1, 6],
	'effect': [doEffect, 7],
	'reading_frame': [doReadingFrame, 8],
	'confidence' :[doConfidence, 9]
}

def filterfunction (data, columnPos): return data[columnPos['confidence']] in ['high', 'medium']

def buildjson(assembly, cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='mutationVector'
	J['dataSubtype'] = 'somatic structural variant'
	J['label'] = 'Arriba fusion'
	J['assembly'] = assembly
	J['cohort'] = cohort
	J['version'] = datetime.date.today().isoformat()
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

filearray_convert.fileArrayFusionToXena(suffix, columns, inputdir, output, columnfunctions, paired_columnfunctions, 
	filterfunction, sample_mapping_file)

buildjson(assembly, cohort, output)