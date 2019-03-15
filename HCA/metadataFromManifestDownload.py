import string, json, sys, os, fnmatch
from sets import Set
sys.path.insert(0, os.path.dirname(__file__))

cell_suspension_file = "cell_suspension_0.json"

bio_metadatafiles =[
	"cell_suspension_0.json",
	"organoid_*.json",
	"cell_line_*.json",
	"specimen_from_organism_*.json",
	"donor_organism_*.json"
]

datafileShortTitle = {
	"cell_suspension": "cell_suspension",
	"cell_line" : "cell_line",
	"organoid" : "organoid",
	"specimen_from_organism": "specimen",
	"donor_organism": "donor"
}


def mtx_output(feature_dic, cell_suspension_id, fout):
	for feature in feature_dic:
                value = map(str, feature_dic[feature])
                value.sort()
		value = string.join(value, '; ')
		fout.write(cell_suspension_id + '\t' + feature + '\t' + value +'\n')

def addToDic (dic, key, feature_data):
	if key not in dic:
		dic[key] =[]
	if feature_data not in dic[key]:
		dic[key].append(feature_data)
	return dic

def parseMeta_to_mtx (indir):
	dic ={}

	# go through bio metadata files
	for file_pattern in bio_metadatafiles:
		files = fnmatch.filter(os.listdir(indir), file_pattern)
		for file in files:
			input = indir + "/" + file
			fin = open(input, 'U')
			J = json.loads(fin.read())
			fin.close()

			for key in J.keys():
				if key in  ["describedBy","schema_type","provenance","total_estimated_cells"]:
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
					dic = addToDic (dic, feature, feature_data)
				
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
					dic = addToDic (dic, feature, feature_data)
					#feature_list = processFeature(fout, feature, feature_data, feature_list)
				
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
					dic = addToDic (dic, feature, feature_data)
					#feature_list = processFeature(fout, feature, feature_data, feature_list)
				
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
						elif k in ["biomaterial_id", "biomaterial_name", "biomaterial_description"]:
							prefix = string.join(string.split(file[:-4], '_')[:-1], '_')
                                                        if k == "biomaterial_id":
                                                                feature = datafileShortTitle[prefix] + "_id"
                                                        elif k == "biomaterial_description":
                                                                feature = datafileShortTitle[prefix] + "_description"
                                                        else:
                                                                continue
							feature_data = value[k]
							dic = addToDic (dic, feature, feature_data)
						else:
							feature = k
							feature_data = value[k]
							dic = addToDic (dic, feature, feature_data)

	return dic

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
	feature_dic = parseMeta_to_mtx (dir)

	allFeatures = allFeatures.union(feature_dic.keys())
	cells.append(cell_suspension_id)

	mtx_output(feature_dic, cell_suspension_id, fout)

fout.close()

cellIndex = {} # key: cell value: cell index
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
