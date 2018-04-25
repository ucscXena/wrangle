import string, sys
import random

def writeShuffledData(items, cellOrder, fout):
	fout.write(items[0])
	for index in cellOrder:
		fout.write('\t' + items[index])
	fout.write('\n')

def process(infile, outfile, maxN):
	fin = open(infile,'r')
	fout = open(outfile, 'w')
        line = fin.readline()[:-1]
	header = string.split(line,'\t')
	cellN = len(header)
	cellOrder = range(1, cellN)
        random.shuffle(cellOrder)
        cellOrder = cellOrder[:maxN]
	writeShuffledData(header, cellOrder, fout)
	while (1):
		line = fin.readline()[:-1]
		if line == '':
			break
		data = string.split(line, '\t')
		writeShuffledData(data, cellOrder, fout)
	fin.close()
	fout.close()

if len(sys.argv[:]) < 3 :
	print "python selectRandomGenomicsData.py infile outfile"
	print
	sys.exit()
	
infile = sys.argv[1]
outfile = sys.argv[2]
if len(sys.argv[:]) == 4:
        maxN = int(sys.argv[3])
process(infile, outfile, maxN)
