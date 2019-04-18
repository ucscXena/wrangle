import string, sys
from sets import Set

gene_exp_file_pattern = ["filtered_gene_bc_matrices_h5.h5", "rsem.genes.results"]

# only check gene expression data in gene_exp_file_pattern
def parseManifest(manifest):
	fin = open(manifest, 'r')
	fin.readline()
	dic ={} # key is bundle uuid value is bundle version
	for line in fin.readlines():
		data = string.split(line,'\t')
		bundle_uuid = data[0]
		bundle_version = data[1]
		filename = data[3]
		# only check gene expression data in gene_exp_file_pattern
		if len(filter(lambda pattern: filename.find(pattern) != -1, gene_exp_file_pattern)) == 0:
			continue
		if bundle_uuid not in dic:
			dic[bundle_uuid]={}
		dic[bundle_uuid]['version'] = bundle_version
		dic[bundle_uuid]['file'] = filename
	fin.close()
	return dic

def checkSameBundles(sameBundles, oldInfo, newInfo):
	needNew = []
	for bundle_uuid in sameBundles:
		if oldInfo[bundle_uuid]['version'] != newInfo[bundle_uuid]['version'] or \
			oldInfo[bundle_uuid]['file'] != newInfo[bundle_uuid]['file']:
			needNew.append(bundle_uuid)
	return needNew

def output(updatemanifest, add, newmanifest):
	fin = open(newmanifest, 'r')
	fout = open(updatemanifest, 'w')
	
	fout.write(fin.readline())

	for line in fin.readlines():
		data = string.split(line,'\t')
		bundle_uuid = data[0]
		if bundle_uuid in add:
			fout.write(line)
	fin.close()
	fout.close()

def outputList (removeBundlesFile, remove):
	fout = open(removeBundlesFile, 'w')
	for bundle_uuid in remove:
		fout.write(bundle_uuid + '\n')
	fout.close()

if len(sys.argv[:]) != 5:
	print "python metadataFromManifestDownload.py oldManifest newManifest addManifest removeBundleList"
	print
	sys.exit()

oldmanifest = sys.argv[1]
newmanifest = sys.argv[2]
addemanifest = sys.argv[3]
removeBundlesFile = sys.argv[4]

oldInfo = parseManifest(oldmanifest)
newInfo = parseManifest(newmanifest)

oldBundles = Set(oldInfo.keys())
newBundles = Set(newInfo.keys())

#new bundle needs to download
add = newBundles.difference(oldBundles) # newBundles - oldBundles

#old bundle needs to remove
remove = oldBundles.difference(newBundles) # oldBundles - newBundles

#same bundle needs update
sameBundles = newBundles.intersection(oldBundles)
needNew = Set(checkSameBundles(sameBundles, oldInfo, newInfo))

add = add.union(needNew)

print "add", len(add)
print add
print "remove", len(remove)
print remove

if len(add) != 0:
	output(addemanifest, add, newmanifest)

if len(remove) != 0:
	outputList (removeBundlesFile, remove)
