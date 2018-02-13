import string, sys, os

DATAfile = "data.tsv"
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

def filter(infile, good_cells, good_genes, output):
	fin = open(infile, 'r')
	fout = open(output, 'w') 
	#cells
	line = fin.readline()
	cells = string.split(line[:-1],'\t')
	good_cols =[]
	fout.write("cell")
	for i in range (1, len(cells)):
		if cells[i] in good_cells:
			good_cols.append(i)
			fout.write("\t" + cells[i])
	fout.write('\n')

	while 1:
		line = fin.readline()
		if line == '':
			break
		data = string.split(line, '\t')
		gene = data[0]
		fout.write(gene)
		if gene not in good_genes:
			continue
		for i in range(1, len(data)):
			if i in good_cols:
				fout.write('\t' + data[i])
		fout.write('\n')

	fin.close()
	fout.close()

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

output = "data_c"+ str(per_cell_gene_count_min)+"_g"+ str(per_gene_cell_count_min)
filter(data_dir + cellStatfile, good_cells, good_genes, output)
