import string, sys
import sets

def allCells (big_genomic):
	fin = open(big_genomic, 'r')
	cells = string.split(fin.readline()[:-1],'\t')[1:]
	fin.close()
	return cells

def fill_small (small_clin, all_cells, big_clin):
	fin = open(small_clin ,'r')

	#header
	line = fin.readline()
	N = string.split(line,'\t')

	#small_clin_cells
	small_clin_cells =[]
	while 1:
		line = fin.readline()
		if line == '':
			break
		cell = line[:line.find('\t')]
		small_clin_cells.append(cell)
	fin.close()

	s = sets.Set(all_cells)
	t = sets.Set(small_clin_cells)
	missing = s.difference_update(t) #return set s after removing elements found in t

	os.system("cp " + small_clin +" " + big_clin)
	fout = open(big_clin, 'a')
	for cell in missing:
		fout.write(cell + '\t'* (N-1) + '\n')
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