import string, sys
import h5py
import numpy as np


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
    new_indptr = np.zeros(counter_indptr_size, dtype=int)

    N = len(indptr) -1 ### ?

    for i in range (0, counter_indptr_size):
        new_data.append( np.empty(N, dtype=int))
        new_indices.append( np.empty(N, dtype=int))

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    for i in range (0, N):
        if i % 500 == 0:
            print "old", i
        indices_range = indices[indptr[i]:indptr[i+1]]
        data_range = data[indptr[i]:indptr[i+1]]
        for index in range (0, len(indices_range)):
            j = indices_range[index]
            value = data_range[index]
            new_data[j][new_indptr[j]] = value
            new_indices[j][new_indptr[j]] = i
            new_indptr[j] += 1

    for i in range (0, len(new_indptr)):
        new_data[i].resize(new_indptr[i])
        new_indices[i].resize(new_indptr[i])

    fout =open(output,'w')
    fout.write("sample\t"+string.join(barcodes,'\t')+'\n')

    for i in range (0, counter_indptr_size):
        gene = genes[i]
        data_range = new_data[i]
        indices_range = new_indices[i]
        values = np.zeros(N, dtype=int)
        for j in range (0, len(indices_range)):
            index = indices_range[j]
            value = data_range[j]
            values[index] = value
        fout.write(gene+'\t'+string.join(map(lambda x: str(x), values),'\t')+'\n')
    fout.close()

# [start:end)
def h5_to_xena (output, data, indices, indptr, counter_indptr_size, genes, barcodes, start, end):
    fout =open(output,'w')
    if start ==0:
        fout.write("sample\t"+string.join(barcodes,'\t')+'\n')

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    for i in range (start, end):
        gene = genes[i]
        indices_range = indices[indptr[i]:indptr[i+1]]
        data_range = data[indptr[i]:indptr[i+1]]
        values = np.zeros(counter_indptr_size, dtype=int)
        for j in range (0, len(indices_range)):
            index = indices_range[j]
            value = data_range[j]
            values[index] = value
        fout.write(gene+'\t'+string.join(map(lambda x: str(x), values),'\t')+'\n')
    fout.close()

if __name__ == "__main__" and len(sys.argv[:])!= 4 and len(sys.argv[:])!= 6:
    print "pyton h5_xena.py h5_input group_name tsv_output start_list(inclusive) end_list(exlusive)\n"
    print "pyton h5_xena.py h5_input group_name tsv_output\n"
    sys.exit()

matrix_h5 = sys.argv[1]
output = sys.argv[3]
group = sys.argv[2]

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

assert(len(indptr) - 1 == colN)
counter_indptr_size = rowN

#optional start and end
if len(sys.argv[:]) == 6:
    start = int(sys.argv[4])
    end = int(sys.argv[5])
    assert (end < len(indptr))
else:
    start = 0
    end = len(indptr) -1 ### total

h5_to_xena (output, data, indices, indptr, counter_indptr_size, genes, barcodes, start, end)

