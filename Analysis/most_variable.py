import sys, os
import pandas as pd
import numpy as np

Na_value = 0 # NA value
withData_threshold = 0.8 # minimum % of cells with data

def var_count(gMatrixfile1, Na_value = 0, withData_threshold = 0.8):
	df1 = pd.read_csv(gMatrixfile1, sep = '\t', index_col = 0)
	total = df1.shape[1]
	total * withData_threshold

	df1.replace(Na_value, np.nan, inplace= True)
	var = df1.var(axis=1, numeric_only= True)
	non_NA = df1.count(axis=1, numeric_only= True)
	r = pd.concat([var, non_NA], axis=1, keys=['var','count'])

	r.dropna(inplace=True)
	r = r[r['count'] >= total * withData_threshold] 
	print (total * withData_threshold)
	r.sort_values(by=['var'], inplace= True)
	print (r)
	return r

if len(sys.argv[:]) != 4:
	print ("python most_variable.py gMatrix zero_value stats_outoput\n")
	print ("All zero_value will be treated as NA, not used to compute variance.\n")
	print ("Only keep features with at least 80% data != NA, sorted by the variance in the output.\n")
	sys.exit(1)

gMatrixfile1 = sys.argv[1]
Na_value = float(sys.argv[2]) # NA value
output = sys.argv[3]

r = var_count (gMatrixfile1, Na_value, withData_threshold)
r.to_csv(output, sep ='\t')