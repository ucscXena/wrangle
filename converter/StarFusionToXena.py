import sys, os, json, datetime
import requests
sys.path.insert(0,os.path.dirname(__file__))
import filearray_convert

suffix = 'STAR.fusion'
columns = ['LeftBreakpoint', 'RightBreakpoint', 
	'LeftGene', 'RightGene',
	'LeftBreakDinuc', 'RightBreakDinuc',
	'annots', 'PROT_FUSION_TYPE', 'FFPM']

compliment = {
	'A' : 'T',
	'T' : 'A',
	'G' : 'C',
	'C' : 'G'
}

def parsePos(genomicPos):
	Chr, pos, strand = genomicPos.split(':')
	if Chr[0:3].upper() != 'CHR':
		Chr = 'chr' + Chr
	return {
		'Chr': Chr,
		'pos': pos,
		'strand': strand
	}

def getAssembly():
	return sys.argv[4]

def doChr1 (data, columnPos): return parsePos(data[columnPos['LeftBreakpoint']])['Chr']
def doPos1 (data, columnPos): return parsePos(data[columnPos['LeftBreakpoint']])['pos']
def doRef1 (data, columnPos):
	strand = parsePos(data[columnPos['LeftBreakpoint']])['strand']
	DNA = data[columnPos['LeftBreakDinuc']][0]	
	DNA = DNA.upper()
	if strand == '-':
		DNA = compliment[DNA]
	return DNA

def doAlt1 (data, columnPos):
	ref = doRef1(data, columnPos)
	LeftBreakpoint = parsePos(data[columnPos['LeftBreakpoint']])
	RightBreakpoint = parsePos(data[columnPos['RightBreakpoint']])
	
	strandLeft = LeftBreakpoint['strand']
	strandRight = RightBreakpoint['strand']
	if strandLeft == '+' and strandRight == '+':
		d1 = 'downstream' # partner 
		d2 = 'upstream' # partner 
	elif strandLeft == '+' and strandRight == '-':
		d1 = 'downstream' # partner
		d2 = 'downstream' # partner
	elif strandLeft == '-' and strandRight == '+':
		d1 = 'upstream' # partner 
		d2 = 'upstream' # partner 
	elif strandLeft == '-' and strandRight == '-':
		d1 = 'upstream' # partner
		d2 = 'downstream' # partner

	bp2 = RightBreakpoint['Chr'] + ':' + RightBreakpoint['pos']

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

def doGene1(data, columnPos): return data[columnPos['LeftGene']].split('^')[0]
def doAltGene1(data, columnPos): return	data[columnPos['RightGene']].split('^')[0]
def doEffect (data, columnPos): 
	r = json.loads(data[columnPos['annots']])
	return ', '.join(r)
def doReadingFrame (data, columnPos): return data[columnPos['PROT_FUSION_TYPE']]
def doFFPM (data, columnPos): return data[columnPos['FFPM']]

def doChr2 (data, columnPos): return parsePos(data[columnPos['RightBreakpoint']])['Chr']
def doPos2 (data, columnPos): return parsePos(data[columnPos['RightBreakpoint']])['pos']

def doRef2 (data, columnPos):
	strand = parsePos(data[columnPos['RightBreakpoint']])['strand']
	DNA = data[columnPos['RightBreakDinuc']][1]
	DNA = DNA.upper()
	if strand == '-':
		DNA = compliment[DNA]
	return DNA

def doAlt2 (data, columnPos):
	ref = doRef2(data, columnPos)
	LeftBreakpoint = parsePos(data[columnPos['LeftBreakpoint']])
	RightBreakpoint = parsePos(data[columnPos['RightBreakpoint']])
	
	strandLeft = LeftBreakpoint['strand']
	strandRight = RightBreakpoint['strand']
	if strandLeft == '+' and strandRight == '+':
		d2 = 'downstream' # partner 
		d1 = 'upstream' # partner 
	elif strandLeft == '+' and strandRight == '-':
		d2 = 'downstream' # partner
		d1 = 'downstream' # partner
	elif strandLeft == '-' and strandRight == '+':
		d2 = 'upstream' # partner 
		d1 = 'upstream' # partner 
	elif strandLeft == '-' and strandRight == '-':
		d2 = 'upstream' # partner
		d1 = 'downstream' # partner

	bp2 = LeftBreakpoint['Chr'] + ':' + LeftBreakpoint['pos']

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

def doGene2(data, columnPos): return data[columnPos['RightGene']].split('^')[0]
def doAltGene2(data, columnPos): return data[columnPos['LeftGene']].split('^')[0]

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
	'FFPM' :[doFFPM, 9]
}

paired_columnfunctions = {
	'chr': [doChr2, 0],
	'start': [doPos2, 1],
	'end': [doPos2, 2],
	'reference': [doRef2, 3],
	'alt': [doAlt2, 4],
	'gene': [doGene2, 5],
	'altGene': [doAltGene2, 6],
	'effect': [doEffect, 7],
	'reading_frame': [doReadingFrame, 8],
	'FFPM' :[doFFPM, 9]
}

def filterfunction (data, columnPos): return float(data[columnPos['FFPM']]) > 0.1

def buildjson(assembly, cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='mutationVector'
	J['dataSubtype'] = 'gene fusion SV (RNA-seq)'
	J['label'] = 'Star-Fusion'
	J['assembly'] = assembly
	J['cohort'] = cohort
	J['version'] = datetime.date.today().isoformat()
	J['wrangling_procedure'] = "Only fusion variants with FFPM > 0.1 is kept"
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
