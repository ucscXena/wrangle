import sys,string,os

non_silent_cutoff =10

typeDic={
    "Nonsense_Mutation":1,
    "stop_gained":1,
    "Frame_Shift_Del":2,
    "Frame_Shift_Ins":2,
    "frameshift_variant":2,
    "Splice_Site":2,
    "splice_acceptor_variant":4,
    "splice_donor_variant":4,
    "splice_region_variant":4,
    "missense_variant":5,
    "Missense_Mutation":5,
    "missense":5,
    "Nonstop_Mutation":6,
    "stop_lost":6,
    "start_lost":6,
    "start_gained":6,
    "initiator_codon_variant":6,
    "De_novo_Start_OutOfFrame":6,
    "Translation_Start_Site":6,
    "De_novo_Start_InFrame":6,
    "5_prime_UTR_premature_start_codon_gain_variant":6,
    "disruptive_inframe_deletion":7,
    "In_Frame_Del":7,
    "inframe_deletion":7,
    "In_Frame_Ins":8, 
    "inframe_insertion":8,
    "Indel":8,
    "RNA":9,
    "non_coding_exon_variant": 9,
    "exon_variant":10, #non-coding 
    "Silent":100 ,##
    "synonymous_variant":100, ##
    "stop_retained_variant":100, ##
    "3_prime_UTR_variant":100, ##
    "5_prime_UTR_variant":100, ##
    "5'Flank":100,##
    "3'Flank":100,##
    "5'UTR":100,##
    "3'UTR":100,##
    "downstream_gene_variant":100, ##
    "upstream_gene_variant":100,##
    "intergenic_region":100, ##
    "intron_variant": 100, ##
    "Intron":100, ##
    "IGR":100 ##
    }



#xena data file
GENECOL = 4
MutTYPECOL = 7
SAMPLECOL = 0
ALLGENE_FILE = "/inside/home/jzhu/scripts/name2/RefGene.name2"

def allGeneList(file):
    fin = open(file,'r')
    allGenes =[]
    for line in fin.readlines():
        line = string.strip(line)
        if line =="":
            continue
        if line[0]=="#":
            continue
        allGenes.append(line)
    fin.close()
    return allGenes

def process (file, samples, genes, dic):
    fin =open(file, 'r')
    fin.readline()
    
    for line in fin.readlines():
        data =string.split(line[:-1],'\t')
        gene = data[GENECOL]
        sample = data[SAMPLECOL]
        if gene =="":
            continue
        try:
            mtype = typeDic[data[MutTYPECOL]]
        except:
            continue
        if sample not in samples:
            samples.append(sample)
        if gene not in genes:
            genes.append(gene)
        if gene not in dic:
            dic[gene]={}
        if sample not in dic[gene]:
            dic[gene][sample]=mtype
        else:
            if mtype < dic[gene][sample]:
                dic[gene][sample]=mtype
    fin.close()

def output (outfile, samples, genes, dic):
    fout=open(outfile,'w')
    fout.write("sample\t"+string.join(samples,"\t")+"\n")
    for gene in dic:
        fout.write(gene)
        for sample in samples:
            try:
                mtype = dic[gene][sample]
                if mtype < non_silent_cutoff:
                    fout.write("\t1")
                else:
                    fout.write("\t0")
            except:
                fout.write("\t0")
        fout.write("\n")
    fout.close()


if len(sys.argv[:]) < 3:
    print "python mafToMatrix.py matrixOUTPUT  mafTypeFilen\n"
    sys.exit()

outfile = sys.argv[1]
files = sys.argv[2:]
samples=[]
genes=[]
dic={}

allGeneList = allGeneList(ALLGENE_FILE)
for gene in allGeneList:
    dic[gene]={}
        
for file in files:
    process (file, samples, allGeneList, dic)

output(outfile, samples, genes, dic)
