#require datamash

import string, sys, os
import math

def DESEQ2_normalization (input, size_factor_output, output):
	fin = open(size_factor_output, 'r')
	fin.readline()
	data = string.split(fin.readline(),'\t')
	fin.close()	
	size_factor = map(lambda x: float(x),data)
	
	colN = len(size_factor) + 1
	
	fin = open(input, 'r')
	fout = open(output, 'w')
	fout.write(fin.readline())
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line[:-1])
		gene = data[0]
		fout.write(gene)
		for i in range(1, colN):
			value = float(data[i])
			if value > 0.0:
				value = value / size_factor[i-1]
			fout.write("\t" + str(value))
		fout.write('\n')
	fin.close()
	fout.close()

def Geometric_Mean_normalization (infile, size_factor_output):
	fin = open(infile, 'r')
	tmpfile = ".out"
	fout = open(tmpfile, 'w')
	
	line = fin.readline()
	data = string.split(line,'\t')
	fout.write(string.join(data[1:],'\t'))
	colN = len(data)
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line[:-1])

		count = 0
		logTotal = 0
		values = []
		for i in range(1, colN):
			value = float(data[i])
			if value > 0.0:
				count = count +1
				logTotal = logTotal + math.log(value)
			elif value < 0:
				print "error: value can not be less than zero"
				sys.exit()
			values.append(value)
		
		if count >=1:
			geometric_mean = math.exp(logTotal / float(count))	
		else:
			geometric_mean = 1

		values = map(lambda x: x/geometric_mean, values)
		for value in values[:-1]:
			if value != 0.0:
				fout.write(str(value)+'\t')
			else:
				fout.write('NA\t')
		if values[-1] != 0.0:
			fout.write(str(values[-1])+'\n')
		else:
			fout.write('NA\n')
	fin.close()
	fout.close()

	os.system("datamash --narm -H median 1-" + str(colN-1) + " < " + tmpfile + " > " + size_factor_output)


if len(sys.argv[:])!=4:
	print "python DESEQ2_normalization.py input_matrix size_factor_output output_matrix\n"
	sys.exit()

input = sys.argv[1]
size_factor_output = sys.argv[2]
output = sys.argv[3]

Geometric_Mean_normalization (input, size_factor_output)
DESEQ2_normalization (input, size_factor_output, output)

