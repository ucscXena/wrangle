
# coding: utf-8

# ## enter sample

# In[ ]:

Nof1_sample = raw_input('Enter sample name (e.g. 09-3-B1): ') or "09-3-B1"
print Nof1_sample


# ## enter gene

# In[ ]:
import re
genes = raw_input('Enter a single or a list of gene names (e.g. PTEN or PTEN TP53): ') or "PTEN,TP53"
genes = re.split( '\W+', genes)
print genes


# ## Run - results at the bottom

# In[ ]:

Nof1_hub = "https://itomic.xenahubs.net"
Nof1_dataset = "latestCCI_EXP_G_TPM_log"

import xena_datasetlist

comparison_list = [
    xena_datasetlist.TCGA_BRCA_tumors,
    xena_datasetlist.TCGA_TNBC
]

import itomic_Nof1

for gene in genes:
    itomic_Nof1.itomic_Nof1(Nof1_sample, gene, Nof1_hub, Nof1_dataset, comparison_list)

itomic_Nof1.itomic_legend()

