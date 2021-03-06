import string, sys, os, json
import xenaPython as xena

if len(sys.argv[:])!= 2:
	print "python probeMapAudit.py xena/files/dir"
	print "genomic_matrix_file_with_probeMap  if_it_is_in_reference_hub"
	print
	sys.exit()

datadir = sys.argv[1]
reference_probemaps =  xena.probemap_list(xena.PUBLIC_HUBS['reference'])

probemapList =[]
for map in reference_probemaps:
	probemapList.append(map['name'])

for dirpath, dirnames, filenames in os.walk(datadir):	# dirpath, dirnames, filenames
	for file in filenames:
		if file [-5:] == ".json":
			jsonFile = dirpath +'/'+ file
			fin = open(jsonFile, 'U')
			meta = json.load(fin)
			probemap = None
			if 'probeMap' in meta:
				probemap = meta['probeMap']
			elif ':probeMap' in meta:
				probemap = meta[':probeMap']
			if probemap and probemap [0] == '/':
				probemap = probemap[1:]
			#print jsonFile, probemap, 
			if probemap:
				print jsonFile, probemap in probemapList
			#else:
			#	print
