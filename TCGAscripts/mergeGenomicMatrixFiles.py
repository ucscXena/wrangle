import string,os,sys

def header (samples, infile):
    fin= open(infile,'r')
    #header, samples
    newSamples = string.split(string.strip(fin.readline()),'\t')[1:]
    for sample in newSamples:
        if sample not in samples:
            samples[sample]= len(samples)
    fin.close()
    return

def process(genes, samples, dataMatrix, infile):
    maxLength= len(samples)

    fin= open(infile,'r')
    #header 
    newSamples = string.split(string.strip(fin.readline()),'\t')
    
    while 1:
        line = fin.readline()[:-1]
        if line =="":
            break
        data = string.split(line,"\t")
        gene = data[0]
        if gene not in genes:
            genes[gene]= len(genes)
            l=[]
            for i in range (0, maxLength):
                l.append("NA")
            dataMatrix.append(l)

        x = genes[gene]
        for i in range (1, len(data)):
            sample = newSamples[i]
            y = samples[sample]
            dataMatrix[x][y]= data[i]

    fin.close()
    return

def outputMatrix(dataMatrix, samples, genes, outfile):
    fout = open(outfile,"w")
    maxLength= len(samples)
    sList=[]
    for i in range (0, maxLength):
        sList.append("")
    for sample in samples:
        pos =samples[sample]
        sList[pos] = sample

    fout.write("sample")
    for sample in sList:
        fout.write("\t"+sample)
    fout.write("\n")

    for gene in genes:
        fout.write(gene)
        for sample in sList:
            value = dataMatrix[genes[gene]][samples[sample]]
            fout.write("\t"+value)
        fout.write("\n")
    fout.close()
    return

if __name__ == '__main__' :
    if len(sys.argv[:]) <4:
        print "python mergeFilesByColumn.py output inputfile(s)"
        print "**********memory intensive, not for very genomic data with hugo number of probes"
        print "this is merging data A+B=C\n"
        sys.exit()

    inFiles = sys.argv[2:]
    print inFiles
    outfile = sys.argv[1]
    print outfile

    genes={}
    samples={}
    dataMatrix=[]

    for infile in inFiles:
        header (samples, infile)

    for infile in inFiles:
        process(genes, samples, dataMatrix, infile)

    outputMatrix(dataMatrix, samples, genes, outfile)
