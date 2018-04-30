# filter out bad gene and bad cells
# Seurat Normalization: a global-scaling normalization method LogNormalize that normalizes the gene expression measurements
# for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and 
# log-transforms ( log1p in R: the natural logarithm of the given value plus one) the result.

# Seurat Scale:  z-scored 

import string, sys, math

DATAfile = "matrix.tsv"
cellStatfile = "cell_statistics"
geneStatfile = "gene_statistics"
filterNormScalefile = "matrix.filter.norm.scale.tsv"
perCellGeneCutoff = 500
perGeneCellCutoff = 5

def read_cell_stats_totalgene (cell_stat_file):
	fin = open(cell_stat_file)
	fin.readline()
	cell_stats={}
	badCells = {}
	while 1:
		line = fin.readline()[:-1]
		if line == '':
			break
		data = string.split(line, '\t')
		cell = data[0]
		cell_total_gene = data[1]
		if cell_total_gene < perCellGeneCutoff:
			badCells[cell]
		else:
			cell_stats[cell] = cell_total_gene
	return cell_stats, badCells

def read_gene_stats_totalcell (geneStatfile):
	fin = open(geneStatfile)
	fin.readline()
	gene_stats={}
	badGenes = {}
	while 1:
		line = fin.readline()[:-1]
		if line == '':
			break
		data = string.split(line, '\t')
		gene = data[0]
		gene_total_cell = data[1]
		if gene_total_cell < perCellGeneCutoff:
			badGenes[cell]
		else:
			gene_stats[gene] = gene_total_cell
	return gene_stats, badGenes

def scaleFunction(x, average, varSQRT):
	if x == 'NA':
		return x
	return (x - average) / varSQRT

def process (infile, outfile, cell_stats, badCells, gene_stats, badGenes):
	fin = open(infile, 'r')
	fout = open(outfile, 'w')
	
	header = fin.readline()
	colN = len(header)
	cell_pos ={}
	bad_pos ={}
	for i in range (1, colN):
		cell = header[i]
		if cell in badCells:
			bad_pos[i] = cell
		else:
			cell_pos[i] = cell

	fout.write(header)
	while 1:
		line = fin.readline()[:-1]
		if line == '':
			break
		data = string.split(line, '\t')
		gene = data[0]

		if gene in badGenes:
			continue

		fout.write(gene)
		total = 0
		count = 0
		x2total = 0

		# seurate normalization
		for i in range (1, colN):
			if i in bad_pos:
				continue
			try:
				value = float(data[i])
				value = math.log(value / cell_stats[cell_pos[i]] * 10 * 1000 +1)
				total = total + value
				x2total = x2total + value * value
				count = count + 1
			except:
				value = 'NA'
				
			data[i]= value

		# seurate scaling
		average = total / float(count)
		var = x2total / float(count) - average * average
		varSQRT = math.sqrt(var)
		#for value in map(lambda x: (x - average) / varSQRT, data[1:]):  #map(lambda p: myFunc(p, additionalArgument), pages)
		for value in map(lambda x: scaleFunction(x, average, varSQRT), data[1:])
			fout.write('\t' + str(value))
		fout.write('\n')

	fin.close()
	fout.close()

if len(sys.argv[:]) != 2:
	print "python seurat_filter_norm_scale.py data_dir outfile"
	print "Expect input data is called matrix.tsv"
	print
	sys.exit()

data_dir = sys.argv[1]
if data_dir[-1] !="/":
    data_dir = data_dir +'/' 

infile = data_dir + DATAfile
cell_stats, badCells = read_cell_stats_totalgene(cellStatfile)
gene_stats, badGenes = read_gene_stats_totalcell (geneStatfile)
outfile = sys.argv[2]

process (infile, outfile, cell_stats, badCells, gene_stats, badGenes)