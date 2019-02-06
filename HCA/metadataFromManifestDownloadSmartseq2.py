import string, json, sys, os
from sets import Set

cell_suspension_file = "cell_suspension_0.json"

datafiles =[
	"cell_suspension_0.json",
	"cell_line_0.json",
	"specimen_from_organism_0.json",
	"donor_organism_0.json"
]

datafileShortTitle = {
	"cell_line_0.json" : "cell_line",
	"cell_suspension_0.json": "cell_suspension",
	"specimen_from_organism_0.json": "specimen",
	"donor_organism_0.json": "donor"
}

def processFeature(fout, feature, feature_data, feature_list):
	if feature not in feature_list:
		feature_list.append(feature)
		fout.write(cell_id + '\t' + feature + '\t' + str(feature_data) +'\n')
	return feature_list


def parseMeta (indir, cell_id, fout):
	feature_list = []
	# go through metadata files
	for file in datafiles:
		input = indir + "/" + file
		if not os.path.exists(input):
			continue

		fin = open(input, 'U')
		J = json.loads(fin.read())
		fin.close()
		for key in J.keys():
			if key in  ["describedBy","schema_type","provenance","publication"]:
				continue
			value = J[key]
			
			#print input
			#print key, value

			if (not type(value) is list) and (not type(value) is dict) :
				'''
				"cell_line_type": "stem cell",
				'''
				feature = key
				feature_data = value
				feature_list = processFeature(fout, feature, feature_data, feature_list)
			elif type(value) is list and type(value[0]) is dict and value[0].has_key("text"):
				'''
				"selected_cell_type": [
		        {
		            "text": "inhibitory interneurons",
		            "ontology": "CL:0000498",
		            "ontology_label": "inhibitory interneuron"
		        }
	    		]
	    		'''
				feature = key
				feature_data = value[0]["text"]
				feature_list = processFeature(fout, feature, feature_data, feature_list)
			
			elif type(value) is dict and value.has_key("text"):
				'''
				    "cell_type": {
			        "text": "embryonic stem cell",
			        "ontology": "CL:0002322",
			        "ontology_label": "embryonic stem cell"
			    }
				'''
				feature = key
				feature_data = value["text"]
				feature_list = processFeature(fout, feature, feature_data, feature_list)
			elif type(value) is dict:
				'''
				"biomaterial_core": {
				    "biomaterial_id": "D26Dp14B02",
				    "biomaterial_name": "hESC-derived inhibitory interneurons",
				    "biomaterial_description": "hESC-derived inhibitory interneurons at day 26",
				    "ncbi_taxon_id": [
				        9606
				    ],
				    "genotype": "DCX-Citrine"
				}
				'''
				for k in value:
					if type (value[k]) is list:
						continue
					if k == "biomaterial_id":
						feature = datafileShortTitle[file] + '_' + k
					else:
						feature = k
					feature_data = value[k]
					feature_list = processFeature(fout, feature, feature_data, feature_list)
	return feature_list

if len(sys.argv[:]) != 4:
	print "python metadataFromManifestDownloadSmartseq2.py inputdir mtx_type_output clinicalMatrix_output"
	print
	sys.exit()

inputdir = sys.argv[1]
mtx_type_output = sys.argv[2]
clinicalMatrix_output = sys.argv[3]

allFeatures = Set([])
cells = Set([])

fout = open(mtx_type_output, 'w')
for subdir in os.listdir(inputdir):
	if not os.path.isdir(subdir):
		continue
	dir = inputdir + "/" + subdir

	# get the id for the gene expression file, for smart-seq2 currently, that's the id used in the matrix file
	fin = open(dir + "/" + cell_suspension_file, 'U')
	J = json.loads(fin.read())
	fin.close()
	expressionfile = J["provenance"]["document_id"]
	cell_id = expressionfile

	feature_list = parseMeta (dir, cell_id, fout)

	allFeatures = allFeatures.union(feature_list)
	cells = cells.union([cell_id])

fout.close()

cellIndex ={} # key: cell value: cell index
cells = list(cells)
for i in range(0, len(cells)):
	cellIndex[cells[i]] =  i

featureIndex ={} # key: feature value: feature index
allFeatures = list(allFeatures)
for i in range(0, len(allFeatures)):
	featureIndex[allFeatures[i]] =  i

matrix =[]
for i in range (0,len(cells)):
	matrix.append([])
	for j in range (0, len(allFeatures)):
		matrix[i].append('')


fin = open(mtx_type_output, 'U')
while 1:
	line = fin.readline()
	if line == '':
		fin.close()
		break
	cell, feature, value = string.split(line[:-1],'\t')
	matrix[cellIndex[cell]][featureIndex[feature]] = value

fout = open(clinicalMatrix_output, 'w')
fout.write('sample' + '\t' + string.join(allFeatures, '\t') + '\n')
for cell in cells:
	fout.write(cell)
	for feature in allFeatures:
		fout.write('\t' + matrix[cellIndex[cell]][featureIndex[feature]])
	fout.write('\n')
fout.close()