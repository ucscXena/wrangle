import sys, os
sys.path.insert(0,os.path.dirname(__file__))

import filearray_convert

suffix = 'abundance.tsv$'
columns = ["est_counts", "tpm"]
maxNum = 500

if len(sys.argv[:]) < 3:
	print ("python KallistoToXena.py inputFileDir output_prefix sample_id_mapping(optional,id file_id)")
	print ()
	sys.exit()

inputdir = sys.argv[1]
output_prefix = sys.argv[2]
if len(sys.argv[:]) > 3:
	sample_mapping_file = sys.argv[3]
else:
	sample_mapping_file = None

filearray_convert.fileArrayToXena(suffix, columns, maxNum, inputdir, output_prefix, sample_mapping_file)