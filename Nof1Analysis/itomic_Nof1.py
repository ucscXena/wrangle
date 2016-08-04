import string, sys
from Nof1_functions import *

#itomic specific
def get_itomic_Data (Nof1_sample, gene, hub, dataset):
    print 'Sample: ', Nof1_sample
    print 'Gene:', gene

    itomic_samples = dataset_samples(hub, dataset)
    if Nof1_sample not in itomic_samples:
        print (Nof1_sample, "not in dataset latestCCI_EXP_G_TPM_log")
        print ("See samples:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host=https%3A%2F%2Fitomic.xenahubs.net%3A443&dataset=latestCCI_EXP_G_TPM_log&label=gene%20expression%20TPM&allSamples=true")
        print
        sys.exit()

    itomic_Data = {}
    values= Probe_values(hub, dataset, itomic_samples, gene)

    for i in range(0,len(itomic_samples)):
        itomic_Data[itomic_samples[i]]= values[i] #log(TPM+0.001)

    Nof1_value = itomic_Data[Nof1_sample]
    Nof1_TPM = math.pow(2, Nof1_value)-0.001

    print "log2(TPM):", Nof1_value, "TPM:", '{:.2f}'.format(Nof1_TPM)
    return itomic_Data

def itomic_Nof1(Nof1_sample, gene, Nof1_hub, Nof1_dataset, comparison_list):
    itomic_Data = get_itomic_Data (Nof1_sample, gene, Nof1_hub, Nof1_dataset)
    Nof1_value = itomic_Data[Nof1_sample]

    for item in comparison_list:
        hub = item["hub"]
        dataset = item["dataset"]
        samples = item["samples"]
        name = item["name"]

        values = Gene_values(hub, dataset, samples, gene)
        rank, percentage =  rank_and_percentage (Nof1_value, values)
        print
        print name +" ( n=", len(values), "):"
        print "rank:", rank
        print "percentile:",'{:.2f}%'.format(percentage)
        print map(lambda x: '{:.2f}%'.format(rank_and_percentage(x, values)[1]), itomic_Data.values())
