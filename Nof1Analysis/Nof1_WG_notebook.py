# In[ ]:

import sys, string

import Nof1_functions
Nof1_item = {
    "hub" : "https://itomic.xenahubs.net",
    "dataset" : "latestCCI_EXP_G_TPM_log",
    "mode" : "probe",
    "name" : "itomic_Nof1",
    "label" : "itomic_Nof1",
    "log2Theta" : 0.001,
    "unit": "log(TPM)"
}


# In[ ]:

Nof1_sample = raw_input('Enter sample name (e.g. 10-3-B1 or ALL): ') or "ALL"

if Nof1_sample == "ALL":
    Nof1_item["samples"] = Nof1_functions.dataset_samples( Nof1_item["hub"], Nof1_item["dataset"])
else:
    if (Nof1_functions.checkSamples (Nof1_sample, Nof1_item["hub"], Nof1_item["dataset"])):
        sys.exit()
    else:
        Nof1_item["samples"]= [Nof1_sample]

print Nof1_item["samples"]


# # check samples  -- no need 

# # enter gene 
genes = raw_input('Enter a single or a list of gene names (e.g. PTEN, AR, or ALL): ') or "ALL"

if Nof1_item["mode"] == "probe":
    if genes =="ALL":
        genes = Nof1_functions.dataset_fields( Nof1_item["hub"],  Nof1_item["dataset"])
    else:
        genes = string.split(genes,',')
    genaname_mapping ={}

# In[ ]:
#genes = filter(lambda x: x!='', re.split(';|,| |\n', genes))
#new_genes =[]
#new_genes = [string.strip (genes[0])]
#for i in range (1, len(genes)):
#    gene = string.strip(genes[i])
#    if gene[0] =="(" and gene[-1] ==")":
#        new_genes[-1] = new_genes[-1] + " (" + string.strip(gene[1:-1]) +")"
#    else:
#        new_genes.append(gene)
#genes = new_genes
if len(genes) >10:
    print "genes:", genes[:10],"..."
else:
    print "genes:", genes[:10]

# # gene name mapping

# In[ ]:

#genaname_mapping ={
#    "CTLA-4" : "CTLA4",
#    "LAG-3" : "LAG3",
#    "LIV-1" : "SLC39A6",
#    "PD-L1" : "CD274",
#    "PDL1" : "CD274",
#    "PD-L2" : "PDCD1LG2",
#    "PDL2" : "PDCD1LG2",
#    "TROP2" : "TACSTD2",
#    "p16" : "CDKN2A",
#    "p18" : "CDKN2C",
#    "MLL2" :"KMT2D",
#    "CD105" : "ENG",
#    "YB1" : "YBX1"
#}


# # check gene name

# In[ ]:

#def cleanGeneName_Funtion (originalLable):
#    return string.strip(string.split(originalLable,'(')[0])

#if (Nof1_functions.checkFields(genes, genaname_mapping, Nof1_item["hub"], Nof1_item["dataset"], cleanGeneName_Funtion)):
#    sys.exit()
#else:
#    print "pass"


# # Enter output file name

# In[ ]:

print "samples:",Nof1_sample
outputfile = raw_input('Enter output file name (e.g. ' + Nof1_sample +")") or Nof1_sample
outputfile = "Results_Folder/" + outputfile

# ## Run

# In[ ]:


# In[ ]:
import xena_datasetlist

comparison_list = [
    xena_datasetlist.TCGA_TNBC,
    xena_datasetlist.TCGA_BRCA_tumors,
    #xena_datasetlist.TCGA_Breast_Basal,
    #xena_datasetlist.TCGA_Breast_Her2,
    #xena_datasetlist.TCGA_Breast_LumA,
    #xena_datasetlist.TCGA_Breast_LumB,
    #xena_datasetlist.TCGA_Breast_Adjacent_Normal,
]

import itomic_WG_Nof1
for comparison_item in comparison_list:
    itomic_WG_Nof1.itomic_Nof1(Nof1_item, genes, genaname_mapping, comparison_item, outputfile + "_vs_"+ comparison_item["fileLabel"]) 

