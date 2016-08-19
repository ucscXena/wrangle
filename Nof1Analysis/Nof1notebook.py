
# coding: utf-8

# ## enter sample

# In[ ]:

Nof1_sample = raw_input('Enter sample name (e.g. 09-3-B1): ') or "09-3-B1"
print Nof1_sample


# In[ ]:

import sys

import Nof1_functions
Nof1_item = {
    "hub" : "https://itomic.xenahubs.net",
    "dataset" : "latestCCI_EXP_G_TPM_log",
    "mode" : "probe",
    "name" : "itomic_Nof1",
    "label" : "itomic_Nof1",
    "samples" : [Nof1_sample],
    "log2Theta" : 0.001
}


# # check sample

# In[ ]:

if (Nof1_functions.checkSamples (Nof1_sample, Nof1_item["hub"], Nof1_item["dataset"])):
    sys.exit()
else:
    print "pass"


# # enter gene 

# In[ ]:

import re
import string
genes = raw_input('Enter a single or a list of gene names (e.g. PTEN or PTEN,TP53 or a column of gene names copied from a spreadsheet): ') or "PTEN,TP53"
genes = filter(lambda x: x!='', re.split(';|,| |\n', genes))
new_genes =[]
new_genes = [string.strip (genes[0])]
for i in range (1, len(genes)):
    gene = string.strip(genes[i])
    if gene[0] =="(" and gene[-1] ==")":
        new_genes[-1] = new_genes[-1] + " (" + string.strip(gene[1:-1]) +")"
    else:
        new_genes.append(gene)
genes = new_genes
print genes


# # gene name mapping

# In[ ]:

genaname_mapping ={
    "CTLA-4" : "CTLA4",
    "LAG-3" : "LAG3",
    "LIV-1" : "SLC39A6",
    "PD-L1" : "CD274",
    "PDL1" : "CD274",
    "PD-L2" : "PDCD1LG2",
    "PDL2" : "PDCD1LG2",
    "TROP2" : "TACSTD2",
    "p16" : "CDKN2A",
    "p18" : "CDKN2C"
}


# # check gene name

# In[ ]:

def cleanGeneName_Funtion (originalLable):
    return string.strip(string.split(originalLable,'(')[0])

if (Nof1_functions.checkFields(genes, genaname_mapping, Nof1_item["hub"], Nof1_item["dataset"], cleanGeneName_Funtion)):
    sys.exit()
else:
    print "pass"


# # Enter output file name

# In[ ]:

outputfile = raw_input('Enter output file name (e.g. result.txt): ') or "result.txt"


# ## Run - results at the bottom

# In[ ]:

print genes
import xena_datasetlist

comparison_list = [
    xena_datasetlist.TCGA_TNBC,
    xena_datasetlist.TCGA_BRCA_tumors,
]

import itomic_Nof1
itomic_Nof1.itomic_Nof1(Nof1_item, genes, genaname_mapping, comparison_list, outputfile)

itomic_Nof1.itomic_legend()

