import string, sys, os

def outputMatrix(outfile, matrix, row, column, Garray, Sarray):
    fout = open(outfile, 'w')
    #header
    values=[]
    for j in range (0, column):
        values.append(Sarray[j])
    line ="sample\t"+string.join(values,"\t")+"\n"
    fout.write(line)
    for i in range(0,row):
        values=[]
        gene = Garray[i]
        for j in range (0, column):
            values.append(str(matrix[i][j]))
        line =gene +"\t"+string.join(values,"\t")+"\n"
        fout.write(line)
    fout.close()


def buildMatrix(infile, outfile):
    samples={}
    genes ={}
    matrix =[]

    sampleFile ="tmpSamples"
    geneFile ="tmpGenes"
    os.system("cut -f 1 "+ infile +" | sort |uniq > "+ sampleFile)
    os.system("cut -f 2 "+ infile +" | sort |uniq > "+ geneFile)

    fin = open(sampleFile,'r')
    for line in fin.readlines():
        sample = string.strip(line)
        if sample not in samples:
            p=len(samples)
            samples[sample]=p
    fin.close()

    fin = open(geneFile,'r')
    for line in fin.readlines():
        gene = string.strip(line)
        if gene not in genes:
            p=len(genes)
            genes[gene]=p
    fin.close()

    #gene by sample matrix
    #per gene
    gene_array = []
    for i in range(0,len(samples)):
        gene_array.append("")
        
    for i in range(0,len(genes)):
        matrix.append(list(gene_array))

    fin =open(infile,'r')
    while 1:
        line = fin.readline()[:-1]
        if line =="":
            break
        sample, gene, value = string.split(line,'\t')
        value = float(value)
        sampleP = samples[sample]
        geneP = genes[gene]
        matrix[geneP][sampleP]= value
    
    Sarray = []
    for i in range(0,len(samples)):
        Sarray.append("")
    for sample in samples:
        Sarray[samples[sample]]=sample

    Garray = []
    for i in range(0,len(genes)):
        Garray.append("")
    for gene in genes:
        Garray[genes[gene]]=gene

    outputMatrix(outfile, matrix, len(genes),len(samples), Garray, Sarray)

if len (sys.argv[:])!=3:
    print "python exp_seq_US.py infile outfile"
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]
buildMatrix(infile, outfile)
