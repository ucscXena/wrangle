import h5py
import scipy.sparse as sp_sparse
import tables
import string, sys, os

def print_attrs(name, obj):
    print name #, obj

#info of the h5 tree obj
def info_h5_datasets (filename):
    hF = h5py.File(filename)
    hF.visititems(print_attrs)

#size of individual ob
def len_obj_h5 (filename, obj):
    hF = h5py.File(filename)
    print len(hF[obj])

def get_h5_info (h5_file):
    info_h5_datasets (h5_file)
    obj = dataset + "/data"
    len_obj_h5(matrix_h5, obj)

def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'], dsets['indices'], dsets['indptr']), shape=dsets['shape'])
            return matrix
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")

#output to mtx format -- still need to add header lines
def output_to_mtx (output, indices, data, start_col, end_col):
    fout = open(output,'w')

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    for col in range (0, colN):
        row_indices = indices[indptr[col]:indptr[col+1]]
        data_indices = data[indptr[col]:indptr[col+1]]
        for i in range (0, len(row_indices)):
            row = row_indices[i]
            value = data_indices[i]
            fout.write(str(row)+'\t'+str(col)+'\t'+str(value)+'\n')
    fout.close()

#output to dense matrix (xena - transposed orientation) format
def output_to_xena_T (output, genes, rowN, row, data, start_col, end_col):
    fout = open(output,'w')
    if start_col == 0:
        fout.write("sample\t"+string.join(genes,'\t')+'\n')

    #is the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.
    l=[]
    for i in range (0, rowN):
        l.append(0)

    for col in range (start_col, end_col):
        sample = barcodes[col]
        row_indices = indices[indptr[col]:indptr[col+1]]
        data_indices = data[indptr[col]:indptr[col+1]]
        fout.write(sample)
        values = l[:]
        for i in range (0, len(row_indices)):
            row = row_indices[i]
            value = data_indices[i]
            values[row] = value
        values = map(lambda x: fout.write('\t'+str(x)), values)
        fout.write('\n')
        print col

    fout.close()

def output_genes (output, genes, gene_names):
    fout = open(output,'w')
    for i in range (0, len(genes)):
        id = genes[i]
        gene_name = gene_names[i]
        fout.write(id+'\t'+ gene_name+'\n')
    fout.close()

def output_barcodes (output, barcodes):
    fout = open(output,'w')
    for i in range (0, len(barcodes)):
        id = barcodes[i]
        fout.write(id+'\n')
    fout.close()

if __name__ == "__main__" and len(sys.argv[:])!=6:
    print "pyton h5_mtx.py h5_input mtx_outputdir dataset_name start_column end_column"
    sys.exit()

matrix_h5 = sys.argv[1]
#matrix_h5 = "1M_neurons_filtered_gene_bc_matrices_h5.h5"
outputdir = sys.argv[2]
dataset = sys.argv[3]
start_col = int(sys.argv[4])
end_col = int(sys.argv[5])
get_h5_info (matrix_h5)

hF = h5py.File(matrix_h5)
indptr = hF[dataset +"/indptr"]
indices = hF[dataset + "/indices"]
data = hF[dataset + "/data"]
genes = hF[dataset + "/genes"]
gene_names = hF[dataset + "/gene_names"]
barcodes = hF[dataset + "/barcodes"]
shape = hF[dataset + "/shape"]
rowN = shape[0]
colN = shape[1]
print rowN, colN

matrix = get_matrix_from_h5(matrix_h5, dataset)

#output to dense matrix (xena - transposed orientation) format
os.system("mkdir " + outputdir)
output_genes (outputdir +"/genes.tsv", genes, gene_names)
output_barcodes (outputdir +"/barcodes.tsv", barcodes)
output_to_mtx (outputdir+"/matrix.mtx", indices, data, start_col, end_col)

