import string, sys, math
import xena_query as xena
import samplelist

def Gene_values(hub, dataset, samples, gene):
    values = xena.dataset_gene_values (hub, dataset, samples, [gene])
    return values[0]["scores"][0]

def Probe_values(hub, dataset, samples, probe):
    values = xena.dataset_probe_values (hub, dataset, samples, [probe])
    return values[0]

def dataset_samples(hub,dataset):
    return xena.dataset_samples(hub, dataset)

def rank_and_percentage(v, v_list):
    v = float(v)
    l = map(float, v_list)
    l.sort()
    l.reverse() #large to small
    K = len(l)
    for i in range(0, K):
        if v > l[i]:
            return i, int(float(i)/K *100)
    return K, "100%"

def Samples_ITOMIC():
    hub = "https://itomic.xenahubs.net"
    dataset = "latestCCI_EXP_G_TPM_log"
    return dataset_samples(hub, dataset)

def Probe_ITOMIC(sample,  gene):
    hub = "https://itomic.xenahubs.net"
    dataset = "latestCCI_EXP_G_TPM_log"
    samples =[sample]
    values = Probe_values(hub, dataset, samples, gene)
    return values

def Gene_TCGA_tumors (gene):
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"

    #custom sample list: current TCGA_pancan_gene_exp_RSEM_tpm data from tumors (non-normals)
    samples = samplelist.TCGA_tumors
    values = Gene_values(hub, dataset, samples, gene)
    return values

def Gene_TCGA_breast_tumors (gene):
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"

    #custom sample list: current TCGA_pancan_gene_exp_RSEM_tpm data from breast tumors (non-normals)
    samples = samplelist.TCGA_BRCA_tumors
    values = Gene_values(hub, dataset, samples, gene)
    return values

def Gene_TCGA_TNBC (gene):
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"

    #custom sample list: current TCGA_pancan_gene_exp_RSEM_tpm data from breast tumors (non-normals)
    samples = samplelist.TCGA_TNBC
    values = Gene_values(hub, dataset, samples, gene)
    return values


if __name__ == "__main__" and len(sys.argv[:])!=3:
    print "python Nof1.py sample gene"
    print
    sys.exit()

Nof1_sample = sys.argv[1]
gene = sys.argv[2]
sampleData={}

samples = Samples_ITOMIC()
if Nof1_sample not in samples:
    print Nof1_sample, "not in dataset latestCCI_EXP_G_TPM_log"
    print "See samples:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host=https%3A%2F%2Fitomic.xenahubs.net%3A443&dataset=latestCCI_EXP_G_TPM_log&label=gene%20expression%20TPM&allSamples=true"
    print
    sys.exit()

for sample in samples:
    Nof1_value = Probe_ITOMIC(sample, gene)[0] #log(TPM+0.001)
    sampleData[sample]= Nof1_value

Nof1_value = sampleData[Nof1_sample]
Nof1_TPM = math.pow(2, Nof1_value)-0.001
print sample
print gene, "log2(TPM):", Nof1_value, "TPM:", '{:.2f}'.format(Nof1_TPM)

values = Gene_TCGA_tumors (gene)
rank, percentage =  rank_and_percentage (Nof1_value, values)
print "TCGA_tumors (n=", len(values), "):", rank, '{:d}%'.format(percentage), map(
    lambda x: '{:d}%'.format(rank_and_percentage(x, values)[1]), sampleData.values())

values = Gene_TCGA_breast_tumors (gene)
rank, percentage =  rank_and_percentage (Nof1_value, values)
print "TCGA_BRCA_tumors (n=", len(values), "):", rank, '{:d}%'.format(percentage), map(
    lambda x: '{:d}%'.format(rank_and_percentage(x, values)[1]), sampleData.values())


values = Gene_TCGA_TNBC (gene)
rank, percentage =  rank_and_percentage (Nof1_value, values)
print "TCGA_TNBC (n=", len(values), "):", rank, '{:d}%'.format(percentage), map(
    lambda x: '{:d}%'.format(rank_and_percentage(x, values)[1]), sampleData.values())
