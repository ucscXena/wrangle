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

def convertExp(expfile, expXenafile, mapping):
	# first line header change
	fin = open(expfile, 'r')
	fout = open(expXenafile, 'w')
	data = fin.readline()[:-1].split('\t')
	fout.write(data[0])
	for i in range (1, len(data)):
		cellkey = data[i]
		intKey = mapping[cellkey]
		fout.write('\t' + intKey)
	fout.write('\n')
	fin.close()
	fout.close()
	# copy rest of the file
	os.system('tail -n +2 ' + expfile + ' >> ' + expXenafile)

def convertCell(cellfile, cellXenafile, mapping):
	fin = open(cellfile, 'r')
	fout = open(cellXenafile, 'w')

	#header lien
	fout.write('xena_cell_id\t')
	fout.write(fin.readline())

	# rest
	while 1:
		line = fin.readline()
		if line == '':
			break
		cellkey = line.split('\t')[0]
		intKey = mapping[cellkey]
		fout.write(intKey+'\t')
		fout.write(line) 
	fin.close()
	fout.close()


if len(sys.argv[:])!= 2:
    print ("python tsv.int.py dataDir")
    sys.exit()

dir = sys.argv[1]

mappingfile= dir + '/mapping'
mapping = parse(mappingfile)

expfile = dir + '/expression.tsv'
expXenafile = dir + '/expression.xena.tsv'
convertExp(expfile, expXenafile, mapping)

cellfile = dir + '/cells.tsv'
cellXenafile = dir + '/cells.xena.tsv'
convertCell(cellfile, cellXenafile, mapping)