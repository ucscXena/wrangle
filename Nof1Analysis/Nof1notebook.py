
# coding: utf-8

# ## enter sample

# In[ ]:

Nof1_sample = raw_input('Enter sample name (e.g. 09-3-B1): ')
print Nof1_sample


# ## enter gene 

# In[ ]:

gene = raw_input('Enter gene name (e.g. PTEN): ')
print gene


# ## Run - results at the bottom

# In[ ]:

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

import itomic_Nof1

def legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "Percentile (%) is the percentile of samples has higher expression than", Nof1_sample+"."
    
itomic_Nof1.itomic_Nof1(Nof1_sample, gene, Nof1_hub, Nof1_dataset, comparison_list)
legend()

