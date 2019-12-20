import os, re, uuid, math

def findColumns(file, columns):
	columnPos ={}
	fin = open(file, 'r')
	# header
	while 1:
		line = fin.readline()
		if line == '':
			break
		if line[0] == '#':
			continue
		data = line.strip().split('\t')
		for i in range (0, len(data)):
			for column in columns:
				if data[i] == column:
					columnPos[column] = i
					break
		break
	fin.close()
	return columnPos

def findVcfColumns(file, columns):
	columnPos ={}
	fin = open(file, 'r')
	# header
	while 1:
		line = fin.readline()
		if line == '':
			break
		if line[0:2] == '##':
			continue
		elif line[0] == '#':
			data = line[1:].strip().split('\t')
			for i in range (0, len(data)):
				for column in columns:
					if data[i] == column:
						columnPos[column] = i
						break
		break
	fin.close()
	return columnPos

def setupDir():
	tmpdir = str(uuid.uuid4())
	os.system('mkdir ' + tmpdir)
	return tmpdir

def collectDir(tmpdir, dataFileList, columns):
	for column in columns:
		os.system("paste " + tmpdir + "/*" + column + ' > ' + tmpdir + '_' + column + ".txt")
		dataFileList[column].append(tmpdir + '_' + column + ".txt")
	os.system("rm -rf " + tmpdir)

def Log2plusTheta(inputFile, outputFile, theta):
    fin = open(inputFile,'r')
    fout = open(outputFile,'w')

    #header
    line = fin.readline()
    fout.write(line)

    #nCOL
    nCOL = len(line.strip().split('\t'))

    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = line[:-1].split('\t')
        if nCOL != len(data):
            print ("Wrong nCOL", data[0])
            sys.exit(1)

        fout.write(data[0])

        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                fout.write("\tNA")
            else:
                value = math.log((float(data[i]) + theta),2)
                value = "%.4f" % (value)
                fout.write("\t"+value)
        fout.write("\n")
    fin.close()
    fout.close()

def finalizeDir(dataFileList, output_prefix, columns, pseudocounts):
	for i in range (0, len(columns)):
		column = columns[i]
		pseudocount = pseudocounts[i]
		files = dataFileList[column]
		output = output_prefix + '_' + column + '.txt'
		tmp_output = output_prefix + '_' + column + '.tmp'

		os.system("paste id " + ' '.join(files) + ' > ' + tmp_output)
		Log2plusTheta(tmp_output, output, pseudocount)
		os.system('rm -f ' + ' '.join(files))
		os.system('rm -f ' + tmp_output)

	os.system('rm -f id')

def sampleInfo(sample_mapping_file):
	dic ={}
	fin = open(sample_mapping_file, 'r')
	for line in fin.readlines():
		id, file = line.strip().split('\t')
		dic[file] = id
		if re.search('.gz$', file) != -1:
			file = re.sub('.gz$', '', file)
			dic[file] = id
	fin.close()
	return dic

def parseTable(input, sample, columnPos, columnfunctions, header, fout):
	fin = open(input, 'r')
	readData =0
	while 1:
		line = fin.readline()
		if line == '':
			break
		if line[0]=='#':
			continue
		if not readData:
			readData = 1
			continue
		data = line[:-1].split('\t')
		fout.write(sample)
		for column in header:
			colfunction = columnfunctions[column][0]
			value = colfunction(data, columnPos)
			fout.write('\t' + value)
		fout.write('\n')
	fin.close()

def parseVcf(input, sample, columnPos, columnfunctions, filterfunction, header, fout):
	fin = open(input, 'r')
	while 1:
		line = fin.readline()
		if line == '':
			break
		if line[0]=='#':
			continue
		data = line[:-1].split('\t')
		if filterfunction(data, columnPos):
			fout.write(sample)
			for column in header:
				colfunction = columnfunctions[column][0]
				value = colfunction(data, columnPos)
				fout.write('\t' + value)
			fout.write('\n')
	fin.close()

def fileArrayExpToXena(suffix, columns, maxNum, inputdir, output_prefix, pseudocounts, sample_mapping_file = None):
	if sample_mapping_file:
		file_sample_info = sampleInfo(sample_mapping_file)
	else:
		file_sample_info = None

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue
		firstfile = re.sub("/$", "", inputdir) + '/'+ file
		break

	columnPos = findColumns(firstfile, columns)
	os.system('cut -f 1 ' + firstfile + " > id")

	dataFileList = {}
	for column in columnPos:
		dataFileList[column] =[]
	tmpdir = setupDir()
	counter = 0

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue
		counter = counter + 1
		filename = re.sub("/$", "", inputdir) + '/'+ file
		for column in columns:
			pos = columnPos[column]
			output = tmpdir + "/" + str(counter) + column
			if file_sample_info:
				sample = file_sample_info[file]
			else:
				print file
				sample = file.split('.')[0]
			os.system("echo " + sample + " > " + output)
			os.system('cut -f ' + str(pos + 1) + ' ' + filename + " | tail -n +2 >> " + output)
		if counter == maxNum:
			collectDir(tmpdir, dataFileList, columns)
			tmpdir = setupDir()
			counter = 0
	collectDir(tmpdir, dataFileList, columns)

	finalizeDir(dataFileList, output_prefix, columns, pseudocounts)

def buildTableHeader(columnfunctions):
	columnHeader = range(0,len(columnfunctions.keys()))
	index = columnfunctions.keys()
	for i in range(0, len(index)):
		column = index[i]
		pos = columnfunctions[column][1]
		columnHeader[pos] = column
	return columnHeader

def fileArrayTableToXena(suffix, columns, inputdir, output, columnfunctions, sample_mapping_file = None):
	if sample_mapping_file:
		file_sample_info = sampleInfo(sample_mapping_file)
	else:
		file_sample_info = None

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue
		firstfile = re.sub("/$", "", inputdir) + '/'+ file
		break

	columnPos = findColumns(firstfile, columns)
	fout = open(output, 'w')
	header = buildTableHeader(columnfunctions)
	fout.write('sample\t' + '\t'.join(header) + '\n')

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue

		filename = re.sub("/$", "", inputdir) + '/'+ file

		if file_sample_info:
			sample = file_sample_info[file]
		else:
			print file
			sample = file.split('.')[0]

		parseTable(filename, sample, columnPos, columnfunctions, header, fout)
	fout.close()

def fileArrayVcfToXena(suffix, columns, inputdir, output, columnfunctions, filterfunction, sample_mapping_file = None):
	if sample_mapping_file:
		file_sample_info = sampleInfo(sample_mapping_file)
	else:
		file_sample_info = None

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue
		firstfile = re.sub("/$", "", inputdir) + '/'+ file
		break

	columnPos = findVcfColumns(firstfile, columns)
	fout = open(output, 'w')
	header = buildTableHeader(columnfunctions)
	fout.write('sample\t' + '\t'.join(header) + '\n')

	for file in os.listdir(inputdir):
		if re.search(suffix, file) == -1:
			continue

		filename = re.sub("/$", "", inputdir) + '/'+ file

		if file_sample_info:
			sample = file_sample_info[file]
		else:
			sample = file.split('.')[0]

		parseVcf(filename, sample, columnPos, columnfunctions, filterfunction, header, fout)
	fout.close()


