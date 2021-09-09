n_components =  3
dens_lambda = 1
columns = ["UMAP_1","UMAP_2","UMAP_3"]

import sys

if len(sys.argv[:]) !=  4:
	print("python pcaTo3DUmap.py  pca_input(no_header) meta_input(first_column_id_same_order_with_header) umap_output\n")
	sys.exit()

import numpy as np, pandas as pd
import umap


input = sys.argv[1]
meta_input = sys.argv[2]
output_file = sys.argv[3]

data = np.loadtxt(input)
meta = pd.read_csv(meta_input, delimiter='\t', index_col=0)

u = umap.UMAP(densmap=True, dens_lambda= dens_lambda,  n_components = n_components).fit(data)
output = pd.DataFrame(u.embedding_, columns = columns, index = meta.index)
output.to_csv(output_file, sep= '\t')
