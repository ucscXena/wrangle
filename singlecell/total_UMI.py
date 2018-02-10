import string, sys, os

def column_N(input):
	fin = open(input,'r')
	N = len(string.split(fin.readline(),'\t'))
        fin.close()
	return N

def clinMatrixOut(input, output):
        fin = open(input,'r')
        fout = open(output,'w')
        fout.write("cell\ttotal_UMI_counts\n")
        cells = string.split(fin.readline()[:-1],'\t')
        values = string.split(fin.readline()[:-1],'\t')
        for i in range (0, len(cells)):
                cell = cells[i][4:-1]
                value = values[i]
                fout.write(cell +'\t'+ value+'\n')
        fin.close()
        fout.close()

if len(sys.argv[:])!=3:
	print "python total_UMI.py xena_matrix_input output"
	sys.exit()

input = sys.argv[1]
output = sys.argv[2]

colN = column_N(input)
os.system("datamash -H sum 2-" + str(colN) + " < "+ input + "> .out")

clinMatrixOut(".out", output)

