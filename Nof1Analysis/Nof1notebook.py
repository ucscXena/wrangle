
# coding: utf-8

# ## enter sample

# In[1]:

Nof1_sample = "10-3-B1"


# ## enter gene 

# In[2]:

gene = "CDKN2A"


# ## Run - results at the bottom

# In[3]:

Nof1_hub = "https://itomic.xenahubs.net"
Nof1_dataset = "latestCCI_EXP_G_TPM_log"

import samplelist

comparison_list = [
    { 
    "hub": "https://toil.xenahubs.net",
    "dataset" : "tcga_RSEM_gene_tpm",
    "samples": samplelist.TCGA_tumors,
    "name": "TCGA_tumors"
    },
    { 
    "hub": "https://toil.xenahubs.net",
    "dataset" : "tcga_RSEM_gene_tpm",
    "samples": samplelist.TCGA_BRCA_tumors,
    "name": "TCGA_BRCA_tumors"
    },
    { 
    "hub": "https://toil.xenahubs.net",
    "dataset" : "tcga_RSEM_gene_tpm",
    "samples": samplelist.TCGA_TNBC,
    "name": "TCGA_TNBC"
    },
    { 
    "hub": "https://toil.xenahubs.net",
    "dataset" : "gtex_RSEM_gene_tpm",
    "samples": samplelist.GTEX_breast,
    "name": "GTEX_breast"
    }]


import sys, getopt
import itomic_Nof1

def input(argv, Nof1_sample, gene):
    try:
        opts, args = getopt.getopt(argv,"hs:g:")
        for opt, arg in opts:
            if opt == '-h':
                print ("python Nof1.py -s <sample> -g <gene>")
                sys.exit()
            if opt in ("-s"):
                Nof1_sample = arg
            if opt in ("-g"):
                gene = arg

    except getopt.GetoptError:
        pass
    return Nof1_sample, gene

def legend():
    print "\nExpression values are sorted from high to low, low rank and percentile mean high expression."
    
if __name__ == "__main__":
    Nof1_sample, gene = input(sys.argv[1:], Nof1_sample, gene)
    itomic_Nof1.itomic_Nof1(Nof1_sample, gene, Nof1_hub, Nof1_dataset, comparison_list)
    legend()

