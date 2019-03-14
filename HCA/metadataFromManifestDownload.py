import string, json, sys, os, os.path
from sets import Set
sys.path.insert(0, os.path.dirname(__file__))

cell_suspension_file = "cell_suspension_0.json"

bio_metadatafiles =[
	"cell_suspension_0.json",
	"specimen_from_organism_0.json",
	"donor_organism_0.json"
]

datafileShortTitle = {
	"cell_line_0.json" : "cell_line",
	"specimen_from_organism_0.json": "specimen",
	"donor_organism_0.json": "donor"
}

def processFeature(fout, feature, feature_data, feature_list):
	if feature not in feature_list:
		feature_list.append(feature)
		fout.write(cell_suspension_id + '\t' + feature + '\t' + str(feature_data) +'\n')
	return feature_list


def parseMeta_to_mtx (indir, cell_suspension_id, fout):
	feature_list = []

	# go through bio metadata files
	for file in bio_metadatafiles:
		input = indir + "/" + file
		
		if not os.path.exists(input):
			continue

		fin = open(input, 'U')
		J = json.loads(fin.read())
		fin.close()

		for key in J.keys():
			if key in  ["describedBy","schema_type","provenance"]:
				continue

			value = J[key]
			
			#print input
			#print key, value

			if (not type(value) is list) and (not type(value) is dict) :
				'''
				"total_estimated_cells": 1800,
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

				"cell_morphology": {
			        "cell_morphology": "Normal",
			        "percent_cell_viability": 92,
			        "cell_viability_method": "Trypan blue / manual haemacytometer (C-chip)"
			    },
				'''
				for k in value:
					if type (value[k]) is list:
						continue
					elif k in ["biomaterial_id", "biomaterial_name"]:
						continue
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
cells = []

fout = open(mtx_type_output, 'w')
for subdir in os.listdir(inputdir):
	dir = inputdir + "/" + subdir

	if not os.path.isdir(dir):
		continue

	# get the id for the cell suspension
	# for smart-seq2, cell suspension id is the id used in the matrix file
	# h5 10xgenomics file, cell suspension id is used as prefix for barcode
        import cell_suspension_id
	cell_suspension_id = cell_suspension_id.cellSuspensionID (dir + "/" + cell_suspension_file)
	feature_list = parseMeta_to_mtx (dir, cell_suspension_id, fout)

	allFeatures = allFeatures.union(feature_list)
	cells.append(cell_suspension_id)

fout.close()

cellIndex ={} # key: cell value: cell index
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
