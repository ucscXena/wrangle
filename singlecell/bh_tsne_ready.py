import string, sys, os

def bh_tsne_ready(input, output):
	os.system("head -n 1 " + input + " | tr '\t' '\n' | tail -n +2 > bh_tsne_cells")
	os.system("tail -n +2 " + input + " | cut -f 2- > .bh_tsne_T")
	dir = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))
	os.system("python "+ dir +"/support/transpose.py .bh_tsne_T " + output)
	os.system("rm .bh_tsne_T")

if len(sys.argv[:])!= 3:
	print "python bh_tsne_ready.py input_matrix output_bhtsne(transpose no headers)"
	print
	sys.exit()

input = sys.argv[1]
output = sys.argv[2]

bh_tsne_ready(input, output)

print "cells are written to bh_tsne_cells"