import string, sys
import random

def writeShuffledData(items, cellOrder, fout):
	fout.write(item[0])
	for index in cellOrder:
		fout.write('\t' + items[index])
	fout.write('\n')

def process(infile, outfile):
	fin = open(infile,'r')
	fout = open(outfile, 'w')
	header = string.split(fin.readline(),'\t')
	cellN = len(header)
	cellOrder = random.shuffle(1, cellN)
	writeShuffledData(header, cellOrder, fout)
	while (1):
		line = fin.readline()
		if line == '':
			break
		data = string.split(fin.readline(),'\t')
		writeShuffledData(data, cellOrder, fout)
	fin.close()
	fout.close()

if len(sys.argv[:])!=3:
	print "python selectRandomGenomicsData.py. infile outfile"
	print

infile = sys.argv[1]
outfile = sys.argv[2]
	