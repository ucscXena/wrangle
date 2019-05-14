import string, sys, os
import sets

def allCells (big_genomic):
	fin = open(big_genomic, 'r')
	cells = string.split(fin.readline()[:-1],'\t')[1:]
	fin.close()
	return cells

def deleteextra_big (big_clin, all_cells, small_clin):
	fin = open(big_clin ,'r')

	#header
	line = fin.readline()
	N = len (string.split(line,'\t'))

	while 1:
		line = fin.readline()
		if line == '':
			break
		cell = line[:line.find('\t')]
		if cell in all_cells:
			fout.write(line)
	fin.close()
	fout.close()

if len(sys.argv[:]) != 4:
	print "python deleteextra_clinicalMatrix_to_full.py in_big_clinical in_big_genomic out_small_clinical"
	print
	sys.exit()

big_clin = sys.argv[1]
big_genomic = sys.argv[2]
small_clin = sys.argv[3]

cells = allCells (big_genomic)

deleteextra_big (big_clin, cells, small_clin)
