import sys, json, datetime
import re

extra_at_beginning = 1
extra_at_end = 6
interval = 2

def parseSampleName(value):
	return value.replace('Log Ratio','').strip()

def process_CPTACC_prot_header(line, N):
	data = line.split('\t')
	l = [data[0]]
	for i in range (extra_at_beginning, N - extra_at_end, interval):
		l.append(data[i])
	return map(parseSampleName, l)

def process_CPTACC_prot_line(line, N):
	data = line.split('\t')
	l = [data[0]]
	for i in range (extra_at_beginning, N - extra_at_end, interval):
		l.append(data[i])
	return l

def process_CPTACC_prot(inputfile, outputfile):
	fin = open(inputfile, 'r')
	fout = open (outputfile, 'w')
	
	line = fin.readline()
	N = len(line.split('\t'))
	r = process_CPTACC_prot_header(line, N)
	fout.write('\t'.join(r) + '\n')

	while 1:
		line = fin.readline()
		if line == '':
			break
		r = process_CPTACC_prot_line (line, N)
		fout.write('\t'.join(r) + '\n')
	fin.close()
	fout.close()

def buildjson(cohort, output):
	output = output + '.json'
	fout = open(output, 'w')
	J = {}
	J['type'] ='genomicMatrix'
	J['dataSubtype'] = 'proteomics'
	J['cohort'] = cohort
	J['unit'] = 'log2(sample/pooled reference)'
	J['label'] = 'protein abundance (Tandem mass tag)'
	J['version'] = datetime.date.today().isoformat()
	json.dump(J, fout, indent = 4)
	fout.close()

if len(sys.argv[:]) != 4:
	print ("python cptaccProteinToXena.py prot_infile outfile cohort\n")
	sys.exit(1)

inputfile = sys.argv[1]
outputfile = sys.argv[2]
cohort = sys.argv[3]

process_CPTACC_prot(inputfile, outputfile)
buildjson(cohort, outputfile)