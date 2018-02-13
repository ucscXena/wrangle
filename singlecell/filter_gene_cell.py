import string, sys, os
import uuid

DATAfile = "matrix.tsv"
cellStatfile = "cell_statistics"
geneStatfile = "gene_statistics"

def check(infile, min_threshold):
	fin = open(infile, 'r')
	fin.readline()
	good_list =[]
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line, '\t')
		id = data[0]
		count = int(data[1])
		if count >= min_threshold:
			good_list.append(id)
	fin.close()
	return good_list

def filter(infile, good_cells, good_genes, data_dir, output):
	fin = open(infile, 'r')
	tmpDir =data_dir + str(uuid.uuid4())
        os.system("mkdir "+ tmpDir)
	ftmp = open(tmpDir + "/.out", 'w')
	
	#genes
	ftmp.write(fin.readline())
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line, '\t')
		gene = data[0]
		if gene not in good_genes:
			continue
		ftmp.write(line)
	ftmp.close()
	fin.close()

	#cells
	fin = open(infile, 'r')
	line = fin.readline()
	cells = string.split(line[:-1],'\t')
	good_cols =[]
	for i in range (1, len(cells)):
		if cells[i] in good_cells:
			good_cols.append(i)
	fin.close()

	#cut
	count = 0
	os.system("cut -f 1 " + tmpDir + '/.out > ' + tmpDir + "/data_" + str(count))

	K=5000 

	for i in range(0, len(good_cols), K):
		count = count + 1
		col_list =[]
		for j in range (i, i+K):
			if j >= len(good_cols):
				break
			col_list.append(str(good_cols[j]))
		os.system("cut -f " + string.join(col_list,',') + ' ' + tmpDir + '/.out >' + tmpDir + "/data_" + str(count))	
	
	#paste
	files = []
	for i in range(0, count):
		files.append(tmpDir+"/data_" + str(count))
	os.system("paste " + string.join(files,' ') + ' > ' + output)

	#clean up
	os.system("rm -rf " + tmpDir)

if len(sys.argv[:])!=4:
	print "python filter_gene_cell.py data_dir cell(per_cell_gene_count_min) gene(per_gene_cell_count_min)"
	sys.exit()

data_dir = sys.argv[1]
if data_dir[-1] !="/":
    data_dir = data_dir +'/'

per_cell_gene_count_min = int(sys.argv[2])
per_gene_cell_count_min = int(sys.argv[3])

good_cells = check(data_dir + cellStatfile, per_cell_gene_count_min)
good_genes = check(data_dir + geneStatfile, per_gene_cell_count_min)

print "Cells:", len(good_cells)
print "Genes:", len(good_genes)

output = "matrix_c"+ str(per_cell_gene_count_min)+"_g"+ str(per_gene_cell_count_min)
filter(data_dir + DATAfile, good_cells, good_genes, data_dir, output)
