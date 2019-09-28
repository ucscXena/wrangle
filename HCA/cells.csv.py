import sys, os

def parse(mappingfile):
	mapping = {}
	fin = open(mappingfile, 'r')
	while 1:
		line = fin.readline()
		if line == '':
			break
		cell_int, cellkey = line.split()
		mapping[cellkey] = cell_int
	fin.close()
	return mapping

def convertCell(cellfile, cellXenafile, mapping):
	fin = open(cellfile, 'r')
	fout = open(cellXenafile, 'w')

	#header line
	fout.write('xena_cell_id\t')
	fout.write(fin.readline().replace(',', '\t'))

	# rest
	while 1:
		line = fin.readline()
		if line == '':
			break
		cellkey = line.split('\t')[0]
		intKey = mapping[cellkey]
		fout.write(intKey+'\t')
		fout.write(line.replace(',', '\t')) 
	fin.close()
	fout.close()

if len(sys.argv[:])!= 2:
    print "python cells.csv.py dataDir"
    sys.exit()

dir = sys.argv[1]

infile = dir + '/cells.csv' 
outfile = dir + '/cells.tsv'
mappingfile= dir + '/mapping'

mapping = parse(mappingfile)
convertCell(infile, outfile, mapping)

# os.system("cat " + infile + " | tr , '\t' > " + outfile)
