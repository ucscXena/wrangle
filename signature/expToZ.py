import sys, string, math



def convertToGeneLevelZscore(expdata, output):
	fin = open(expdata, 'U')
	fout = open(output,'w')

	line = fin.readline()
	fout.write(line)
	N = len(string.split(line[:-1],'\t'))
	count = 0
	while 1:
		count = count +1
		line = fin.readline()
		if line == '':
			break
		data = string.split(line[:-1],'\t')
		gene = data[0]
		# relative expression value. Ei,j = log2(TPMi,j/10+1)
		values = map(lambda x: math.log(float(x)/10 +1, 2), data[1:])
		total = reduce(lambda x, y: x+y, values)
		average = total / (N-1)

		t = 0.0
		for i in range (0, N-1):
			t = (values[i] - average) * (values[i] - average) + t
		sigma = math.sqrt(t/ (N-2))
		#print count, total, average, sigma
		#print gene, average, sigma

		if sigma ==0.0:
			pass
		else:
			Z = map(lambda x: (x-average)/sigma, values) # z score
		fout.write(gene+'\t')
		fout.write(string.join(map(lambda x: str(x), Z),'\t'))
		fout.write('\n')
	fin.close()
	fout.close()
	
if len(sys.argv[:])!= 3:
	print "python expToZ.py input output"
	print
	sys.exit()

expdata = sys.argv[1]
zdata = sys.argv[2]


convertToGeneLevelZscore(expdata, zdata)