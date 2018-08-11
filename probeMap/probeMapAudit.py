import string, sys, os, json
import xenaPython as xena

if len(sys.argv[:])!= 2:
	print "python probeMapAudit.py xena/files/dir"
	print
	sys.exit()

datadir = sys.argv[1]
reference_probemaps =  xena.probemap_list(xena.PUBLIC_HUBS['reference'])

probemapList =[]
for map in reference_probemaps:
	print map['name']
	probemapList.append(map['name'])

for dirpath, dirnames, filenames in os.walk(datadir):	# dirpath, dirnames, filenames
	for file in filenames:
		if file [-5:] == ".json":
			jsonFile = dirpath +'/'+ file
			fin = open(jsonFile, 'U')
			meta = json.load(fin)
			probemap = None
			if 'probemap' in meta:
				probemap = meta['probemap']
			elif ':probemap' in meta:
				probemap = meta[':probemap']
			if probemap and probemap [0] == '/':
				probemap = probemap[1:]
			print jsonFile, probemap, 
			if probemap:
				print probemap in probemapList
			else:
				print
