import string, sys, os
import sets
import gzip, mimetypes

def allCells (big_genomic):
	type = mimetypes.guess_type(big_genomic)[1]
	if type == 'gzip':
		fin = gzip.open(big_genomic, 'rb')
	else:
		fin = open(big_genomic, 'r')

	cells = string.split(fin.readline()[:-1],'\t')[1:]
	fin.close()
	return cells

def fill_small (small_clin, all_cells, big_clin, fill_value):
	fin = open(small_clin ,'r')

	#header
	line = fin.readline()
	N = len (string.split(line,'\t'))

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
	s.difference_update(t) #return set s after removing elements found in t
	missing = s

	os.system("cp " + small_clin +" " + big_clin)
	fout = open(big_clin, 'a')
        fill_value = ('\t'+fill_value) * (N-1)
	for cell in missing:
		fout.write(cell + fill_value + '\n')
	fout.close()

if len(sys.argv[:]) not in [4,5]:
	print "python fill_clinicalMatrix_to_full.py in_small_clinical in_big_genomic(gzip or uncompressed) out_big_clinical fill_in_value(optional)"
	print
	sys.exit()

small_clin = sys.argv[1]
big_genomic = sys.argv[2]
big_clin = sys.argv[3]
fill_value = ''

if len(sys.argv[:]) == 5:
	fill_value = sys.argv[4]


cells = allCells (big_genomic)
fill_small (small_clin, cells, big_clin, fill_value)
