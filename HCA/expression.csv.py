import string, sys, os
import csv
from itertools import izip
import uuid

def parse_gene(genes_file):
    probeID_gene_mapping = {}
    fin = open(genes_file, 'r')
    fin.readline()
    for line in fin.readlines():
        data = line[:-1].split(',') # HCA genes csv
        probe = data[0]
        gene = data[1]
        probeID_gene_mapping[probe] = gene
    fin.close()
    return probeID_gene_mapping

def parse(mappingfile):
    mapping = {}
    fin = open(mappingfile, 'r')
    while 1:
        line = fin.readline()
        if line == '':
            break
        cell_int, cellkey = line.split()
        mapping[cellkey] = cell_int
    fin.close()
    return mapping

def process (infile, outfile, input_delimiter, k, total_cell, probeID_gene_mapping):
    fout=open(outfile,'w')
    fout.close()
    fout=open(outfile,'a')

    # read header count gene number
    fin=open(infile,'rU')
    line = fin.readline()
    totalN = len(string.split(line, input_delimiter))
    fin.close()

    count = 0
    tmpdir = str(uuid.uuid4())
    os.system("mkdir " + tmpdir)

    #first line in output file because mapping is generated using the exp file, it is just sequential nubmers
    fout.write('xena_cell_id')
    for i in range(1, total_cell + 1):
        fout.write('\t' + str(i))
    fout.write('\n')

    # rest columns
    for i in range (1, totalN, k):
        count = count + 1
        print (count)
        tmpinfile = tmpdir + "/" + str(count)
        start = i + 1 #linux cut
        command = "cut -d , -f " + str(start) + "-" + str(start+k-1) +" " + infile + " > " + tmpinfile
        os.system(command)

        fin=open(tmpinfile,'rU')
        a = izip(*csv.reader(fin, delimiter = input_delimiter))
        fin.close()

        #switch gene ids in a
        for row in a:
            probeID = row[0]
            try:
                gene = probeID_gene_mapping[probeID]
            except:
                print (probeID)
                gene = probeID
            fout.write(gene+'\t')
            csv.writer(fout, delimiter="\t", lineterminator="\n").writerow(row[1:])
        os.system("rm -rf " + tmpinfile)

    fout.close()
    os.system("rm -rf " + tmpdir)

# transpose
# switch ensemble to gene
if len(sys.argv[:])!= 3:
    print ("python expression.csv.py dataDir chunk_size(e.g.200 larger the file, smaller the size)")
    sys.exit()

dir = sys.argv[1]
infile = dir + '/expression.csv' 
genes_file = dir + '/genes.csv'
outfile= dir + '/expression.tsv'
mappingfile= dir + '/mapping'

mapping = parse(mappingfile)
total_cell = len(mapping)
k = int(sys.argv[2])
print (total_cell)

probeID_gene_mapping = parse_gene(genes_file)
input_delimiter = ',' # HCA matrix csv
process(infile, outfile, input_delimiter, k, total_cell, probeID_gene_mapping)

