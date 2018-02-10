import string, sys, os

def column_N(input):
	fin = open(input,'r')
	N = len(string.split(fin.readline(),'\t'))
	return N

if len(sys.argv[:])!=3:
	print "python total_cellcolumn.py xena_matrix_input output"
	sys.exit()

input = sys.argv[1]
colN = column_N(input)
for i in range(1,colN):
	os.system("cut -f "+ str(i) + " | datamash sum 1 > .out")