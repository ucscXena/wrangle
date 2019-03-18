import string, sys
import h5py

gencode_gene_probemap = "/mnt/singlecell/xena/files/gencode.v27.annotation.gene.Comprehensive.CHR.probemap"
group = "GRCh38"

if len(sys.argv[:]) != 3:
	print "python HCA_h5_probeMap.py input_example_h5_file output_probeMap"
	print "use gencode.v27.annotation.gene.Comprehensive.CHR.probemap as the base"
	print
	sys.exit()

matrix_h5 = sys.argv[1]

hF = h5py.File(matrix_h5, 'r')
genes = hF[group + "/genes"]
hF.close()