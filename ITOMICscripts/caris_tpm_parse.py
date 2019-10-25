import sys

def process (input, output, sample):
	fin = open(input, 'r')
	fout = open(output, 'w')
	fout.write('Gene\t' + sample +'\n')
	while 1:
		line = fin.readline()
		if line == '':
			break
		if line[:2] =="##":
			continue
		fout.write(line)
		fout.write(fin.read())
		break
	fin.close()
	fout.close()

if len(sys.argv[:]) != 4:
    print ("python caris_tpm_parse.py input output sampleLabel")
    sys.exit()

input = sys.argv[1]
output = sys.argv[2]
sample = sys.argv[3]

process(input,  output, sample)