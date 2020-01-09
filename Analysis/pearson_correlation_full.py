import sys, os
import pandas as pd
import numpy as np

def printBasicInfo(df1):
	print (df1.shape)
	print (df1.columns)
	print (df1.index)
	print (df1.dtypes)

# to find profile that are very very similar pairs, used in matching the same sample across dataset,
# although very close, the profile across different dataset are not identical due to processing difference
def find_very_close_pair(corr): 
	N = corr.shape[0]
	for col in corr.columns:
		r = corr[col].sort_values()
		print (col, r.index[N-2], r[N-2]) # N-1 is the corr =1

def compute_N_by_N_pearson_correlation(gMatrixfile1, Na_value = 0):
	df1 = pd.read_csv(gMatrixfile1, sep = '\t', index_col = 0)

	df1.replace(Na_value, np.nan, inplace= True)

	corr = df1.corr(method = 'pearson')

	return corr



if len(sys.argv[:]) != 4:
	print ("python pearson_correlation_full.py gMatrix zero_value stats_outoput\n")
	print ("All zero_value will be treated as NA, not used to compute correlation.\n")
	sys.exit(1)

gMatrixfile1 = sys.argv[1]
Na_value = float(sys.argv[2]) # NA value
output = sys.argv[3]

corr = compute_N_by_N_pearson_correlation(gMatrixfile1, Na_value)
corr.to_csv(output, sep ='\t')
find_very_close_pair(corr)