import string, json, sys, os

if len(sys.argv[:]) != 2:
	print "python metadataFromManifestDownloadSmartseq2.py inputdir"
	print
	sys.exit()

inputdir = sys.argv[1]

