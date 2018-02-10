import string, sys, os

def column_N(input):
	fin = open(input,'r')
	N = len(string.split(fin.readline(),'\t'))
	return N

if len(sys.argv[:])!=3:
	print "python total_UMI.py xena_matrix_input output"
	sys.exit()

input = sys.argv[1]
output = sys.argv[2]
fout = open(output,'w')
fout.close()

colN = column_N(input)

os.system("datamash -H sum 2-" + str(colN) + " < "+ input + ">>" + output)

