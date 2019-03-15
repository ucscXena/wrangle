import string, sys, os
import h5py
from sets import Set
sys.path.insert(0, os.path.dirname(__file__))

cell_suspension_file = "cell_suspension_0.json"
h5_file = "filtered_gene_bc_matrices_h5.h5"
h5Group = "GRCh38"

def output(barcodes, bundle_uuid, fout):
	for barcode in barcodes:
		fout.write(barcode + '\t' + bundle_uuid +'\n')


if len(sys.argv[:]) != 3:
	print "python barcode_cellsuspension_mapping.py inputdir barcode_cellsuspension_mapping_output"
	print
	sys.exit()

inputdir = sys.argv[1]
mapping_output = sys.argv[2]

fout = open(mapping_output, 'w')
for subdir in os.listdir(inputdir):
	dir = inputdir + "/" + subdir

	if not os.path.isdir(dir):
		continue

	# import cell_suspension_id
	# cell_suspension_id = cell_suspension_id.cellSuspensionID (dir + "/" + cell_suspension_file)
	bundle_uuid = subdir # HCA data has duplicate cell suspension id weird, so use subdir/bundle_uuid for now

	hF = h5py.File(dir + "/" + h5_file)
	barcodes = hF[h5Group + "/barcodes"]
	barcodes = map(lambda x:  bundle_uuid + '_' + x, barcodes)
	output(barcodes, bundle_uuid, fout)

fout.close()
