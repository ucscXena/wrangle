import string, sys, os

DATAfile = "matrix.tsv"
cellStatfile = "cell_statistics"
geneStatfile = "gene_statistics"

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

def writeout(data_dir, cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI):
    fcellStat = open(data_dir + cellStatfile,'w')
    fcellStat.write("cell\ttotal_gene_expressed\ttotal_UMI_counts\n")
    fgeneStat = open(data_dir + geneStatfile,'w')
    fgeneStat.write("gene\ttotal_cell\ttotal_UMI\n")

    for i in range (1, len(cells)):
        cell = cells[i]
        fcellStat.write(cell +'\t'+ str(total_gene[i])+ '\t'+ str(total_UMI[i])+'\n')

    for gene in gene_stats_count:
        fgeneStat.write(gene + '\t' + str(gene_stats_count[gene])+ '\t' + str(gene_stats_UMI[gene]) + '\n')

    fcellStat.close()
    fgeneStat.close()

    os.system("sort -gk 2 " + data_dir + cellStatfile + ' > .out')
    os.system("mv .out " + data_dir + cellStatfile)
    os.system("sort -gk 2 " + data_dir + geneStatfile + ' > .out')
    os.system("mv .out " + data_dir + geneStatfile)

if len(sys.argv[:])!=2:
	print "python matrix_stats.py data_dir"
	sys.exit()

data_dir = sys.argv[1]
if data_dir[-1] !="/":
    data_dir = data_dir +'/' 

colN = column_N(data_dir + DATAfile)
cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI = process (data_dir + DATAfile, colN)

writeout(data_dir, cells, total_gene, total_UMI, gene_stats_count, gene_stats_UMI)
