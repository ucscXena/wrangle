
# coding: utf-8

# ## enter sample

# In[ ]:

Nof1_sample = raw_input('Enter sample name (e.g. 09-3-B1): ') or "09-3-B1"
print Nof1_sample


# In[ ]:

import sys
import Nof1_functions
Nof1_hub = "https://itomic.xenahubs.net"
Nof1_dataset = "latestCCI_EXP_G_TPM_log"


# # check sample

# In[ ]:

if (Nof1_functions.checkSamples (Nof1_sample, Nof1_hub, Nof1_dataset)):
    sys.exit()
else:
    print "pass"


# ## enter gene

genaname_mapping ={
}
# In[ ]:

import re
genes = raw_input('Enter a single or a list of gene names (e.g. PTEN or PTEN TP53): ') or "PTEN,TP53"
genes = filter(lambda x: x!='', re.split(';|,| |\n', genes))
print genes
if (Nof1_functions.checkFields(genes, genaname_mapping, Nof1_hub, Nof1_dataset)):
    sys.exit()
else:
    print "pass"


# ## Run - results at the bottom

# In[ ]:

import xena_datasetlist

comparison_list = [
    xena_datasetlist.TCGA_BRCA_tumors,
    xena_datasetlist.TCGA_TNBC
]

import itomic_Nof1

for gene in genes:
    itomic_Nof1.itomic_Nof1(Nof1_sample, genaname_mapping[gene], Nof1_hub, Nof1_dataset, comparison_list)

itomic_Nof1.itomic_legend()

