import string, sys
from Nof1_functions import *

#itomic specific
def get_itomic_Data (gene, hub, dataset, samples):
    values= Probe_values(hub, dataset, samples, string.upper(gene))
    itomic_Data =dict(zip(samples, values))
    return itomic_Data

def Nof1_output (Nof1_sample, gene, itomic_Data):
    print
    print 'Sample: ', Nof1_sample
    print 'Gene:', gene

    Nof1_value = itomic_Data[Nof1_sample]
    Nof1_TPM = math.pow(2, Nof1_value)-0.001

    print "log2(TPM):", Nof1_value, "TPM:", '{:.2f}'.format(Nof1_TPM)


def itomic_Nof1(Nof1_sample, gene, Nof1_hub, Nof1_dataset, comparison_list):
    itomic_samples = dataset_samples(Nof1_hub, Nof1_dataset)

    if Nof1_sample not in itomic_samples:
        print (Nof1_sample, "not in dataset latestCCI_EXP_G_TPM_log")
        print ("See samples:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host=https%3A%2F%2Fitomic.xenahubs.net%3A443&dataset=latestCCI_EXP_G_TPM_log&label=gene%20expression%20TPM&allSamples=true")
        print
        sys.exit()

    itomic_Data = get_itomic_Data (gene, Nof1_hub, Nof1_dataset, itomic_samples)
    Nof1_output (Nof1_sample, gene, itomic_Data)
    Nof1_value = itomic_Data[Nof1_sample]

    for item in comparison_list:
        hub = item["hub"]
        dataset = item["dataset"]
        samples = item["samples"]
        name = item["name"]
        gene = string.upper(gene)

        values = Gene_values(hub, dataset, samples, gene)
        h_l_values = clean_then_sort_high_to_low (values)
        rank, percentage =  rank_and_percentage (Nof1_value, h_l_values)
        print
        print name +" ( n=", len(h_l_values), "):"
        print "rank:", rank
        for list in map(lambda x: (str(x[0]),'{:d}'.format(rank_and_percentage(x[1], h_l_values)[0])),
            zip(itomic_Data.keys(), itomic_Data.values())):
            print list

def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "Percentile rank (%) is the percentile of samples has higher expression than sample of interest."
    print
