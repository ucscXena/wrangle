import string, sys

def allCells (big_genomic):
	fin = open(big_genomic, 'r')
	cells = string.split(fin.readline()[:-1],'\t')[1:]
	fin.close()
	return cells

def fill_small (small_clin, cells, big_clin):
	fin = open(small_clin ,'r')
	fout = open(big_clin, 'w')

	#header
	line = fin.readline()
	N = string.split(line,'\t')
	fout.write(line)
	#data
	while 1:
		line = fin.readline()
		if line == '':
			break
		cell = line[:line.find('\t')]
		if cell not in cells:
			fout.write(cell + '\t'* (N-1) + '\n')
		else:
			fout.write(line)

	fin.close()
	fout.close()

if len(sys.argv[:]) != 4:
	print "python fill_clinicalMatrix_to_full.py in_small_clinical in_big_genomic out_big_clinical"
	print
	sys.exit()

small_clin = sys.argv[1]
big_genomic = sys.argv[2]
big_clin = sys.argv[3]

cells = allCells (big_genomic)

fill_small (small_clin, cells, big_clin)