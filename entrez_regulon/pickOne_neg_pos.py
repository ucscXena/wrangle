import string,sys

if len(sys.argv[:]) !=3:
	print "python pickOne_neg_pos.py regulator regulon_hugo_infile"
	print
	sys.exit()

def process(regulator, regulon_hugo_infile):
	fin = open(regulon_hugo_infile, 'U')
	fin.readline()
	foutPos = open(regulator+"_pos",'w')
	foutNeg = open(regulator+"_neg",'w')
	for line in fin.readlines():
		data = string.split(line[:-1],'\t')
		if string.upper(data[0]) != string.upper(regulator):
			continue
		MoA = float(data[2])
		if MoA >0:
			foutPos.write(line)
		else:
			foutNeg.write(line)

	fin.close()
	foutNeg.close()
	foutPos.close()

regulator = sys.argv[1]
regulon_hugo_infile = sys.argv[2]

process(regulator, regulon_hugo_infile)