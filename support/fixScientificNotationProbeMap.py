import string,sys

def process (input, output):
	fin = open(input, 'U')
	fout = open(output, 'w')
	header = fin.readline()
	fout.write(header)

	for line in fin.readlines():
		data = string.split(line,'\t')
		assert(len(data) == 6)
		data[3] = str(int(float(data[3])))
		data[4] = str(int(float(data[4])))
		fout.write(string.join(data,'\t'))
	fin.close()
	fout.close()


if len(sys.argv[:])!= 3:
	print "python fixScientificNotationProbeMap.py probeMapIn probeMapOut"
	print
	sys.exit()

input = sys.argv[1]
output = sys.argv[2]

process (input, output)