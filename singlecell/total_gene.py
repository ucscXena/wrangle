import string, sys, os

def column_N(input):
    fin = open(input,'r')
    N = len(string.split(fin.readline(),'\t'))
    fin.close()
    return N

def process (input, colN):
    fin = open(input,'r')
    line = fin.readline()
    cells = string.split(line[:-1],'\t')

    total =[]
    for i in range (0, colN):
    	total.append(0)

    count = 0
    while 1:
        count = count +1
        if count % 500 == 0:
            print count
     	line = fin.readline()
    	if line =='':
    		break
    	values = string.split(line[:-1],'\t')
    	for i in range (1, colN):
    		value = float(values[i])
    		if value != 0:
    			total[i] = total[i] +1

    fin.close()
    return cells, total

def writeout(cells, values, output):
        fout = open(output,'w')
        fout.write("cell\ttotal_gene_counts\n")
        for i in range (1, len(cells)):
                cell = cells[i]
                value = values[i]
                fout.write(cell +'\t'+ str(value)+'\n')
        fout.close()

if len(sys.argv[:])!=3:
	print "python total_gene.py xena_matrix_input output"
	sys.exit()

input = sys.argv[1]
output = sys.argv[2]

colN = column_N(input)
cells, total_list = process (input, colN)

#.system("cat "+ input + "| sed 's/\t0\b//'g | datamash -H count 2-" + str(colN) + "> .out")

writeout(cells, total_list, output)

