import string, sys
import h5py

gencode_gene_probemap = "/mnt/singlecell/xena/files/gencode.v27.annotation.gene.Comprehensive.CHR.probemap"
group = "GRCh38"

def parse(gencode_gene_probemap):
	fin = open(gencode_gene_probemap, 'r')
	fin.readline()

	probeMap ={}
	for line in fin.readlines():
		data = string.split(line, '\t')
		probe = data[0]
		probeMap[probe] = line
	fin.close()

	return probeMap

if len(sys.argv[:]) != 3:
	print "python HCA_h5_probeMap.py input_example_h5_file output_probeMap"
	print "use gencode.v27.annotation.gene.Comprehensive.CHR.probemap as the base"
	print
	sys.exit()

matrix_h5 = sys.argv[1]

hF = h5py.File(matrix_h5, 'r')
genes = hF[group + "/genes"]
hF.close()

baseProbeMap = parse(gencode_gene_probemap)

fout = open(sys.argv[2], 'w')
fout.write("id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")
for gene in genes:
	if gene in baseProbeMap:
		fout.write(baseProbeMap[gene])
	else:
		print "ERROR," gene, "not found"
fout.close()

