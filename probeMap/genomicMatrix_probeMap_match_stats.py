import string, sys
from sets import Set

if len(sys.argv[:]) != 3:
	print "python genomicMatrix_probeMap_match_stats.py genomicMatrix probeMap"
	print
	sys.exit()

def parse(genomicMatrix_OR_probemap_file):
	fin = open(genomicMatrix_OR_probemap_file, 'r')
	fin.readline()

	identifiers = []
	while 1:
		line = fin.readline()
		if line == '':
			break
		pos = line.find('\t')
		if pos == -1:
			print "ERROR"
			sys.exit()
		probe = line[:pos]
		identifiers.append(probe)
	fin.close()

	return identifiers

genomicMatrix_file = sys.argv[1]
probeMap_file = sys.argv[2]

probeMap_ids = Set(parse(probeMap_file))
gMX_identifiers = Set(parse(genomicMatrix_file))

overlap = len(probeMap_ids.intersection(gMX_identifiers))

print "probeMap id:", len(probeMap_ids)
print "genomicMatrix id:", len(gMX_identifiers)
print "overlap:", overlap

if overlap < len(gMX_identifiers):
	print
	print "missing in probeMap"
	print list(gMX_identifiers.difference(probeMap_ids))[:10]
	print

if overlap < len(probeMap_ids):
	print
	print "extra in probeMap"
	print list(probeMap_ids.difference(gMX_identifiers))[:10]
	print
