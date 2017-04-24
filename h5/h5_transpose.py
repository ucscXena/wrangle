import h5py
import string, sys
import itertools
import numpy as np

def print_attrs(name, obj):
    print name, len(obj)

def get_h5_info (h5_file):
    hF = h5py.File(h5_file)
    print hF.keys()
    hF.visititems(print_attrs)

def transpose_h5 (data, indices, indptr, new_indptr_size):
    new_data= []
    new_indices= []
    new_indptr = np.zeros(new_indptr_size, dtype=int)

    N = len(indptr) -1 ### ?

    for i in range (0, new_indptr_size):
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
            new_indptr[j] = new_indptr[j] + 1

    for i in range (0, new_indptr_size):
        new_data[i] = new_data[i][:new_indptr[i]]
        new_indices[i] = new_indices[i][:new_indptr[i]]

    total = 0
    for i in range (0, new_indptr_size):
        total = total + new_indptr[i]
        new_indptr[i] = total
    new_indptr = np.insert(new_indptr,0,0)

    # due to memory size, build data, indices and empty lists one by one
    return_data = np.zeros(data.size, dtype=int)
    return_indices = np.zeros(data.size, dtype=int)
    for i in range (0, new_indptr_size):
        return_data[new_indptr[i]: new_indptr[i+1]]=new_data[i]
        new_data[i] = np.array(0)
        return_indices[new_indptr[i]: new_indptr[i+1]]=new_indices[i]
        new_indices[i] = np.array(0)

    #return [np.concatenate(new_data, axis=0 ), np.concatenate(new_indices, axis=0 ), new_indptr]
    return [return_data, return_indices, new_indptr]

def output_h5 (output, group, data, indices, indptr, shape, genes, gene_names, barcodes):
    f = h5py.File(output,'w')
    g = f.create_group(group)
    g.create_dataset('data', data= data, compression="gzip")
    g.create_dataset('indptr',data= indptr, compression="gzip")
    g.create_dataset('indices',data= indices, compression="gzip")
    g.create_dataset('genes', data = genes, compression="gzip")
    g.create_dataset('gene_names', data = gene_names, compression="gzip")
    g.create_dataset('barcodes', data = barcodes, compression="gzip")
    g.create_dataset('shape', data = shape)
    g.attrs['shape'] = shape
    f.close()


if __name__ == "__main__" and len(sys.argv[:])!=4:
    print "pyton h5_transpose.py h5_input group_name h5_output"
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

new_indptr_size = rowN
data, indices, indptr = transpose_h5 (data, indices, indptr, new_indptr_size)
new_shape = [colN,rowN]

output_h5 (output, group, data, indices, indptr, new_shape, genes, gene_names, barcodes)





