import h5py
import string, sys
import numpy

def print_attrs(name, obj):
    print name, len(obj)

def get_h5_info (h5_file):
    hF = h5py.File(h5_file)
    print hF.keys()
    hF.visititems(print_attrs)

#output to mtx format -- still need to add header lines
def output_to_mtx (output, colN, indices, data):
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

def h5_T_to_xena (output, data, indices, indptr, counter_indptr_size, genes, barcodes):
    new_data= []
    new_indices= []

    for i in range (0, counter_indptr_size):
        new_data.append([])
        new_indices.append([])

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    N = len(indptr) -1 ### ?
    for i in range (0, N):
        indices_range = indices[indptr[i]:indptr[i+1]]
        data_range = data[indptr[i]:indptr[i+1]]
        for index in range (0, len(indices_range)):
            j = indices_range[index]
            value = data_range[index]
            new_data[j].append(value)
            new_indices[j].append(i)

    fout =open(output,'w')
    fout.write("sample\t"+string.join(barcodes,'\t')+'\n')

    l=[]
    for i in range (0, N):
        l.append(0)

    for i in range (0, counter_indptr_size):
        gene = genes[i]
        data_range = new_data[i]
        indices_range = new_indices[i]
        values = l[:]
        for j in range (0, len(indices_range)):
            index = indices_range[j]
            value = data_range[j]
            values[index] = value
        fout.write(gene+'\t'+string.join(map(lambda x: str(x), values),'\t')+'\n')
    fout.close()


def h5_to_xena (output, data, indices, indptr, counter_indptr_size, genes, barcodes):
    fout =open(output,'w')
    fout.write("sample\t"+string.join(barcodes,'\t')+'\n')

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    N = len(indptr) -1 ### ?

    l=[]
    for i in range (0, counter_indptr_size):
        l.append(0)

    for i in range (0, N):
        gene = genes[i]
        indices_range = indices[indptr[i]:indptr[i+1]]
        data_range = data[indptr[i]:indptr[i+1]]
        values = l[:]
        for j in range (0, len(indices_range)):
            index = indices_range[j]
            value = data_range[j]
            values[index] = value
        fout.write(gene+'\t'+string.join(map(lambda x: str(x), values),'\t')+'\n')
    fout.close()


if __name__ == "__main__" and len(sys.argv[:])!=4:
    print "pyton h5_xena.py h5_input group_name tsv_output"
    sys.exit()

matrix_h5 = sys.argv[1]
output = sys.argv[3]
group = sys.argv[2]
get_h5_info (matrix_h5)

hF = h5py.File(matrix_h5)
indptr = hF[group +"/indptr"]
indices = hF[group + "/indices"]
data = hF[group + "/data"]
genes = hF[group + "/genes"]
gene_names = hF[group + "/gene_names"]
barcodes = hF[group + "/barcodes"]
shape = hF[group + "/shape"]
rowN = shape[0]
colN = shape[1]
print "row", rowN
print "col", colN

assert(len(indptr) -1 == colN)

counter_indptr_size = rowN
h5_to_xena (output, data, indices, indptr, counter_indptr_size, genes, barcodes)





