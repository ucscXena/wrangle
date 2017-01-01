import sys,string,os

non_silent_cutoff = 1

typeDic={
    'Nonsense_Mutation': 4,
    'Nonsense': 4,
    'frameshift_variant': 4,
    'Frameshift': 4,
    'stop_gained': 4,
    'Stop Gained': 4,
    'Frame_Shift_Del': 4,
    'Frame_Shift_Ins': 4,
    'Frameshift Deletion': 4,
    'Frameshift Insertion': 4,

    #splice related
    'splice_acceptor_variant': 3,
    'splice_acceptor_variant&intron_variant': 3,
    'splice_donor_variant': 3,
    'splice_donor_variant&intron_variant': 3,
    'SpliceAcceptorDeletion': 3,
    'SpliceAcceptorSNV': 3,
    'SpliceDonorBlockSubstitution': 3,
    'SpliceDonorDeletion': 3,
    'SpliceDonorSNV': 3,
    'Splice_Site': 3,
    'splice_region_variant': 3,
    'splice_region_variant&intron_variant': 3,

    #modify protein
    'missense': 2,
    'non_coding_exon_variant': 2,
    'missense_variant': 2,
    'Missense Variant': 2,
    'Missense_Mutation': 2,
    'Missense': 2,
    'MultiAAMissense': 2,
    'start_lost': 2,
    'start_gained': 2,
    'De_novo_Start_OutOfFrame': 2,
    'Translation_Start_Site': 2,
    'CdsStartSNV': 2,
    'De_novo_Start_InFrame': 2,
    'stop_lost': 2,
    'Stop Lost': 2,
    'Nonstop_Mutation': 2,
    'initiator_codon_variant': 2,
    '5_prime_UTR_premature_start_codon_gain_variant': 2,
    'disruptive_inframe_deletion': 2,
    'disruptive_inframe_insertion': 2,
    'inframe_deletion': 2,
    'Inframe Deletion': 2,
    'InFrameDeletion': 2,
    'inframe_insertion': 2,
    'Inframe Insertion': 2,
    'InFrameInsertion': 2,
    'In_Frame_Del': 2,
    'In_Frame_Ins': 2,
    'Indel': 2,

    #do not modify protein
    'synonymous_variant': 1,
    'Synonymous Variant': 1,
    'Synonymous': 1,
    'Silent': 1,
    'stop_retained_variant': 1,

    #mutations effect we don't know
    'lincRNA': 0,
    'RNA': 0,
    'exon_variant': 0,
    'upstream_gene_variant': 0,
    'downstream_gene_variant': 0,
    "5'Flank": 0,
    "3'Flank": 0,
    "3'UTR": 0,
    "5'UTR": 0,
    '5_prime_UTR_variant': 0,
    '3_prime_UTR_variant': 0,
    'Complex Substitution': 0,
    'intron_variant': 0,
    'intergenic_region': 0,
    'intron': 0,
    'Intron': 0,
    "IGR":0
}


#xena data file
#ALLGENE_FILE = "/inside/home/jzhu/scripts/name2/RefGene.name2"
#ALLGENE_FILE = "/data/TCGA/tcgaDataOneOff/Genomic/PANCAN/genes"
ALLGENE_FILE = "/inside/home/jzhu/cgDataJing/TCGAscripts/gencode_genes"

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

    #header parser
    line = fin.readline()
    data =string.split(line[:-1],'\t')
    
    SAMPLECOL = 0

    for i in range (0, len(data)):
        if data[i]=="gene":
            GENECOL = i
        if data[i]=="effect":
            MutTYPECOL = i
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
            if mtype > dic[gene][sample]:
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
                if mtype > non_silent_cutoff:
                    fout.write("\t1")
                else:
                    fout.write("\t0")
            except:
                fout.write("\t0")
        fout.write("\n")
    fout.close()


if __name__ == '__main__':
    if len(sys.argv[:]) < 3:
        print "python xenaToMatrix.py matrixOut xenaFile(s)\n"
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
