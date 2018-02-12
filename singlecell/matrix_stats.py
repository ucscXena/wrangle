import string, sys, os

cellStatoutput = "cell_statistics"
geneStatoutput = "gene_statistics"

def column_N(input):
    fin = open(input,'r')
    N = len(string.split(fin.readline(),'\t'))
    fin.close()
    return N

def process (input, colN):
    fin = open(input,'r')
    line = fin.readline()
    cells = string.split(line[:-1],'\t')

    total_gene =[]
    total_UMI =[]
    gene_stats_UMI ={}
    gene_stats_count ={}

    for i in range (0, colN):
    	total_gene.append(0)
        total_UMI.append(0)

    count = 0
    while 1:
        count = count +1
        if count % 500 == 0:
            print count
     	line = fin.readline()
    	if line =='':
    		break
    	values = string.split(line[:-1],'\t')
        gene = values[0]
        if gene in gene_stats_UMI:
            print "error, duplicated", gene
        if gene in gene_stats_count:
            print "error, duplicated", gene
        gene_stats_UMI[gene]=0
        gene_stats_count[gene] = 0

    	for i in range (1, colN):
            value = float(values[i])
            if value != 0:
                total_gene[i] = total_gene[i] +1
                gene_stats_count[gene] = gene_stats_count[gene] + 1

            total_UMI[i] = total_UMI[i] + value
            gene_stats_UMI[gene] = gene_stats_UMI[gene] + value
    fin.close()

    return cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI

def writeout(outputdir, cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI):
    fcellStat = open(outputdir + cellStatoutput,'w')
    fcellStat.write("cell\ttotal_gene_expressed\ttotal_UMI_counts\n")
    fgeneStat = open(outputdir + geneStatoutput,'w')
    fgeneStat.write("gene\ttotal_cell\ttotal_UMI\n")

    for i in range (1, len(cells)):
        cell = cells[i]
        fcellStat.write(cell +'\t'+ str(total_gene[i])+ '\t'+ str(total_UMI[i])+'\n')

    for gene in gene_stats_count:
        fgeneStat.write(gene + '\t' + str(gene_stats_count[gene])+ '\t' + str(gene_stats_UMI[gene]) + '\n')

    fcellStat.close()
    fgeneStat.close()

    os.system("sort -gk 2 " + outputdir + cellStatoutput + ' > .out')
    os.system("mv .out > " + outputdir + cellStatoutput)
    os.system("sort -gk 2 " + outputdir + geneStatoutput + ' > .out')
    os.system("mv .out > " + outputdir + geneStatoutput)

if len(sys.argv[:])!=3:
	print "python matrix_stats.py xena_matrix_input outputdir"
	sys.exit()

input = sys.argv[1]
outputdir = sys.argv[2]
if outputdir[-1] !="/":
    outputdir = outputdir +'/' 

colN = column_N(input)
cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI = process (input, colN)

writeout(outputdir, cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI)
