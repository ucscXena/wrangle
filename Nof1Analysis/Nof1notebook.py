
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


# # enter gene

# In[ ]:

import re
genes = raw_input('Enter a single or a list of gene names (e.g. PTEN or PTEN TP53): ') or "PTEN,TP53"
genes = filter(lambda x: x!='', re.split(';|,| |\n', genes))
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
    "TROP2" : "TACSTD2"
}


# # check gene name

# In[ ]:

if (Nof1_functions.checkFields(genes, genaname_mapping, Nof1_hub, Nof1_dataset)):
    sys.exit()
else:
    print "pass"


# ## Run - results at the bottom

# In[ ]:

import xena_datasetlist

comparison_list = [
    xena_datasetlist.TCGA_TNBC,
    xena_datasetlist.TCGA_BRCA_tumors,
]

import string
outputfile = "result.txt"
fout = open(outputfile,'w')
headerList =["gene"]
header2ndList =[""]
for item in comparison_list:
    headerList.append(item["label"]+ ' (n=' + str(len(item["samples"]))+")")
    headerList.append("Range of previous ITOMIC samples vs. " + item["label"])
    header2ndList.append("expression %")
    header2ndList.append("expression %")
headerList.extend(["",""])
header2ndList.extend(["log2(TPM)","TPM"])

fout.write(string.join(headerList,'\t') +'\n')
fout.write(string.join(header2ndList,'\t') +'\n')

import itomic_Nof1
for gene in genes:
    itomic_Nof1.itomic_Nof1(Nof1_sample, gene, genaname_mapping[gene], Nof1_hub, Nof1_dataset, comparison_list, fout)
fout.close()

itomic_Nof1.itomic_legend()

