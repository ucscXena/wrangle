import sys, os
sys.path.insert(0,os.path.dirname(__file__))

import filearray_convert

suffix = 'abundance.tsv$'
columns = ["est_counts", "tpm"]
maxNum = 500

if len(sys.argv[:]) != 3:
	print ("python KallistoToXena.py inputFileDir output_prefix")
	print ()
	sys.exit()

inputdir = sys.argv[1]
output_prefix = sys.argv[2]

filearray_convert.fileArrayToXena(suffix, columns, maxNum, inputdir, output_prefix)