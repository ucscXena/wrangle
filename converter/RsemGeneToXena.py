import sys, os
sys.path.insert(0,os.path.dirname(__file__))

import filearray_convert

suffix = 'genes.results$'
columns = ["expected_count", "TPM", 'FPKM']
pseudocounts = [1, 0.001, 0.001]
maxNum = 500

if len(sys.argv[:]) < 3:
	print ("python RsemGeneToXena.py inputFileDir output_prefix sample_id_mapping(optional, id file_id) probeMapfilename")
	print ()
	sys.exit(1)

inputdir = sys.argv[1]
output_prefix = sys.argv[2]
if len(sys.argv[:]) > 3:
	sample_mapping_file = sys.argv[3]
else:
	sample_mapping_file = None
if len(sys.argv[:]) >4 :
	probeMapfilename = sys.argv[4]
else:
	probeMapfilename = None

filearray_convert.fileArrayToXena(suffix, columns, maxNum, inputdir, output_prefix, pseudocounts, 
	sample_mapping_file, probeMapfilename)
