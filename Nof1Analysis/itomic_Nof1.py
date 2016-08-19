import string, sys
from Nof1_functions import *

#itomic specific
def get_itomic_Data (gene, hub, dataset, samples):
    values= Probe_values(hub, dataset, samples, gene)
    itomic_Data =dict(zip(samples, values))
    return itomic_Data

def itomic_log2TPM_2_TPM (Nof1_value):
    theta =0.001
    return math.pow(2, Nof1_value)- theta

def Nof1_output (Nof1_sample, original_label, gene, itomic_Data):
    print
    print 'Sample: ', Nof1_sample
    print 'Gene:', original_label
    if original_label!= gene:
        print 'HUGO gene name:', gene

    Nof1_value = itomic_Data[Nof1_sample]
    Nof1_TPM = itomic_log2TPM_2_TPM(Nof1_value)

    print "log2(TPM):", Nof1_value, "TPM:", '{:.2f}'.format(Nof1_TPM)


def itomic_Nof1(Nof1_sample, original_label, gene, Nof1_hub, Nof1_dataset, comparison_list, fout):
    itomic_samples = dataset_samples(Nof1_hub, Nof1_dataset)

    itomic_Data = get_itomic_Data (gene, Nof1_hub, Nof1_dataset, itomic_samples)
    Nof1_output (Nof1_sample, original_label, gene, itomic_Data)
    Nof1_value = itomic_Data[Nof1_sample]

    outputList =[gene]

    print Nof1_value,"here"
    for item in comparison_list:
        hub = item["hub"]
        dataset = item["dataset"]
        samples = item["samples"]
        name = item["name"]
        gene = string.upper(gene)

        values = Gene_values(hub, dataset, samples, gene)
        h_l_values = clean_then_sort_high_to_low (values)

        rank, percentage =  rank_and_percentage (Nof1_value, h_l_values)
        outputList.append('{:.2f}%'.format(percentage))

        r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), itomic_Data.values())
        outputList.append('['+string.join(map(lambda x: '{:.2f}%'.format(x[1]), r_and_p_values),', ')+']')

        print
        print name +" ( n=", len(h_l_values), "):"
        print "rank:", rank
        print map(lambda x: x[0], r_and_p_values)

        print "expression level percentile:", '{:.2f}%'.format(percentage)

        for list in zip(itomic_Data.keys(), r_and_p_values):
            print list[0], list[1][0], '{:.2f}%'.format(list[1][1])

    outputList.append(str(Nof1_value))
    outputList.append(str(itomic_log2TPM_2_TPM(Nof1_value)))
    fout.write(string.join(outputList,'\t') +'\n')


def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "expression level percentile is the percentile of samples has lower expression than sample of interest."
    print "High expression level percentile means higher expression."
    print
