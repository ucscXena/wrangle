import h5py
import string, sys
import numpy

def print_attrs(name, obj):
    print name, len(obj)

def get_h5_info (h5_file):
    hF = h5py.File(h5_file)
    print hF.keys()
    hF.visititems(print_attrs)

def transpose_h5 (data, indices, indptr, new_indptr_size):
    new_data= []
    new_indices= []
    new_indptr = []

    for i in range (0, new_indptr_size):
        new_data.append([])
        new_indices.append([])

    #the standard CSC representation
    #where the row indices for column i are stored in indices[indptr[i]:indptr[i+1]] and
    #their corresponding values are stored in data[indptr[i]:indptr[i+1]].
    #If the shape parameter is not supplied, the matrix dimensions are inferred from the index arrays.

    N = len(indptr) -1 ### ?
    for i in range (0, N):
        print i
        indices_range = indices[indptr[i]:indptr[i+1]]
        data_range = data[indptr[i]:indptr[i+1]]
        for index in range (0, len(indices_range)):
            j = indices_range[index]
            value = data_range[index]
            new_data[j].append(value)
            new_indices[j].append(i)

    data =[]
    indices =[]
    indptr = []
    index = 0
    for i in range (0, new_indptr_size):
        indptr.append(index)
        index = index + len(new_data[i])
        data.extend(new_data[i])
        new_data[i]=0
        indices.extend(new_indices[i])
        new_indices[i]=0
    indptr.append(index)

    return [data, indices, indptr]

def output_h5 (output, group, data, indices, indptr, shape, genes, gene_names, barcodes):
    f = h5py.File(output,'w')
    g = f.create_group(group)
    g.create_dataset('data', data= data, compression="gzip")
    g.create_dataset('indptr',data=indptr, compression="gzip")
    g.create_dataset('indices',data=indices, compression="gzip")
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





