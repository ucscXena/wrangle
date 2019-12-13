import os, re, uuid

def findColumns(file, columns):
	columnPos ={}
	fin = open(file, 'r')
	# header
	data = fin.readline().strip().split('\t')
	for i in range (0, len(data)):
		for column in columns:
			if data[i] == column:
				columnPos[column] = i
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

def finalizeDir(dataFileList, output_prefix, columns):
	for column in columns:
		files = dataFileList[column]
		os.system("paste id " + ' '.join(files) + ' > ' + output_prefix + '_' + column + '.txt')
		os.system('rm -f ' + ' '.join(files))
	os.system('rm -f id')


def fileArrayToXena(suffix, columns, maxNum, inputdir, output_prefix):
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
			sample = file.split('.')[0]
			os.system("echo " + sample + " > " + output)
			os.system('cut -f ' + str(pos + 1) + ' ' + filename + " | tail -n +2 >> " + output)
		if counter == maxNum:
			collectDir(tmpdir, dataFileList, columns)
			tmpdir = setupDir()
			counter = 0
	collectDir(tmpdir, dataFileList, columns)

	finalizeDir(dataFileList, output_prefix, columns)