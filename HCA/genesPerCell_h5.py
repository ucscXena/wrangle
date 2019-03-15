import string, sys, os
import h5py

cell_suspension_file = "cell_suspension_0.json"
h5_file = "filtered_gene_bc_matrices_h5.h5"
h5Group = "GRCh38"

def genesPerCell_h5 (hF):
	group = h5Group

	indptr = hF[group +"/indptr"]
	indices = hF[group + "/indices"]
	data = hF[group + "/data"]
	genes = hF[group + "/genes"]
	gene_names = hF[group + "/gene_names"]
	barcodes = hF[group + "/barcodes"]
	shape = hF[group + "/shape"]
	rowN = shape[0] # number of genes
	colN = shape[1] # number of samples
	counter_indptr_size = rowN
	    
	#the standard CSC representation
	#where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
	#their corresponding values are stored in data[indptr[i]:indptr[i+1]].
	#If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

	sampleIndex = range(0, colN)
	genesInCell_list = map(lambda i: len(data[indptr[i]:indptr[i+1]]), sampleIndex)
	
	return genesInCell_list

def stats(numbers, fout):
	fout.write('genes_per_cell\tcumulative\n')

	N = len(numbers)
	current_value = 0.0
	for i in range (0, N):
		if numbers[i] > current_value:
			cumulative = i / float(N)
			fout.write(str(current_value) + '\t' + str(cumulative) +'\n')
				current_value = numbers[i]
	fout.write(str(current_value) + '\t1.0\n')
	
if len(sys.argv[:]) != 3:
	print "python genesPerCell_h5.py inputdir output"
	print
	sys.exit()

inputdir = sys.argv[1]
output = sys.argv[2]

for subdir in os.listdir(inputdir):
	dir = inputdir + "/" + subdir

	if not os.path.isdir(dir):
		continue

	if not os.path.exists(dir + "/" + cell_suspension_file):
		continue

	numbers = []
	hF = h5py.File(dir + "/" + h5_file)
	genesInCell_list = genesPerCell_h5 (hF)
	numbers.extend(genesInCell_list)

fout = open(output, 'w')
numbers.sort()
stats(numbers, fout)
fout.close()