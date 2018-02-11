import string, sys, os

geneoutput = "gene_counts_per_sample"
UMIoutput = "UMI_counts_per_sample"
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
    gene_stats ={}

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
        if gene in gene_stats:
            print "error, duplicated", gene
        gene_stats[gene]=0

    	for i in range (1, colN):
            value = float(values[i])
            if value != 0:
                total_gene[i] = total[i] +1
            total_UMI[i] = total_UMI[i] + value
            gene_stats[gene] = gene_stats[gene] + value
    fin.close()

    for gene in gene_stats:
        gene_stats[gene] = gene_stats[gene]/colN

    return cells, total_gene, total_UMI, gene_stats

def writeout(outputdir, cells, total_gene, total_UMI, gene_stats):
    fgene = open(outputdir + geneoutput,'w')
    fgene.write("cell\ttotal_gene_counts\n")
    fUMI = open(outputdir + UMIoutput,'w')
    fUMI.write("cell\ttotal_UMI_counts\n")
    fgeneStat = open(outputdir + geneStatoutput,'w')
    fgeneStat.write("gene\taverage_UMI\n")

    for i in range (1, len(cells)):
        cell = cells[i]
        value = total_gene[i]
        fgene.write(cell +'\t'+ str(value)+'\n')
        value = total_UMI[i]
        fUMI.write(cell +'\t'+ str(value)+'\n')

    for gene in gene_stats:
        fgeneStat.write(gene + '\t' + str(gene_stats[gene])+'\n')

    fgene.close()
    fUMI.close()
    fgeneStat.close()

if len(sys.argv[:])!=3:
	print "python matrix_stats.py xena_matrix_input outputdir"
	sys.exit()

input = sys.argv[1]
outputdir = sys.argv[2]
if outputdir[-1] !="/":
    outputdir = outputdir +'/' 

colN = column_N(input)
cells, total_gene, total_UMI, gene_stats = process (input, colN)

writeout(outputdir, cells, total_gene, total_UMI, gene_stats)
