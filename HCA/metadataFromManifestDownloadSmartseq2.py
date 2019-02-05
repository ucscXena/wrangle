import string, json, sys, os

if len(sys.argv[:]) != 2:
	print "python metadataFromManifestDownloadSmartseq2.py inputdir"
	print
	sys.exit()

inputdir = sys.argv[1]

for subdir in os.listdir(inputdir):
	if os.path.isdir(subdir):
		pass
	else:
		print subdir

