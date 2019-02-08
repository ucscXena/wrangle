import string, json, sys, os

cell_suspension_file = "cell_suspension_0.json"

tmpdatadir = "tmp"
rsem_gene_dir = "rsem_gene"
rsem_gene_file_suffix = "_rsem.genes.results"
rsem_gene_Pos = {
	"expected_count": 5,  # 1-based index to use with bash
	"TPM": 6,
	"FPKM": 7
}

def buildTmpDataDir (tmpdatadir):
	os.system("rm -rf " + tmpdatadir)
	os.system("mkdir " + tmpdatadir)
	os.system("mkdir " + tmpdatadir + "/" + rsem_gene_dir)
	for unit in rsem_gene_Pos.keys():
		os.system("mkdir " + tmpdatadir + "/" + rsem_gene_dir +"/" + unit)

def buildFirstColumn(indir, tmpdatadir, cell_id):
	input = indir + "/" + cell_id + rsem_gene_file_suffix
	output = tmpdatadir + "/" + rsem_gene_dir + "/firstColumn"
	os.system("cut -f 1 " + input + " > " + output)

def parseGenomics (indir, tmpdatadir, cell_id):
	input = indir + "/" + cell_id + rsem_gene_file_suffix
	for unit in rsem_gene_Pos.keys():
		pos = rsem_gene_Pos[unit]
		outputfile = tmpdatadir + "/" + rsem_gene_dir + "/" + unit + "/" + cell_id
		os.system("echo '" + cell_id + "' > " + outputfile )
		os.system("cut -f " + str(pos) + " " + input + " | tail -n +2 >> " + outputfile)

def mergeData(tmpdatadir):
	for unit in rsem_gene_Pos.keys():
		firstColumn = tmpdatadir + "/" + rsem_gene_dir + "/firstColumn"

		inputfilesdir = tmpdatadir + "/" + rsem_gene_dir + "/" + unit
		inputfiles = map(lambda x : inputfilesdir + "/" + x, os.listdir(inputfilesdir))
		count = 0
		outputfiles = []
		for i in range (0, len(inputfiles), 500):
			outputfile = rsem_gene_dir + "_" + unit + '_' + str(count)
			outputfiles.append(outputfile)
			os.system("paste " + string.join(inputfiles[i: i+500], ' ')  + " > " + outputfile)
			count = count + 1
		outputfile = rsem_gene_dir + "_" + unit
		os.system("paste " + firstColumn + " " + string.join(outputfiles, ' ') + " > " + outputfile)
		os.system("rm " + string.join(outputfiles, ' '))

if len(sys.argv[:]) != 2:
	print "python compileData.py inputdir"
	print
	sys.exit()

inputdir = sys.argv[1]

buildTmpDataDir(tmpdatadir)

count = 0
for subdir in os.listdir(inputdir):
	if not os.path.isdir(subdir):
		continue
	if subdir == tmpdatadir:
		continue
	dir = inputdir + "/" + subdir

	count = count +1
	print count, dir
	# get the id for the gene expression file, for smart-seq2 currently, that's the id used in the matrix file
	fin = open(dir + "/" + cell_suspension_file, 'U')
	J = json.loads(fin.read())
	fin.close()
	expressionfile = J["provenance"]["document_id"]
	cell_id = expressionfile
	
	if count == 1:
		buildFirstColumn(dir, tmpdatadir, cell_id)

	parseGenomics (dir, tmpdatadir, cell_id)

mergeData(tmpdatadir)