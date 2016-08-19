import string, sys
from Nof1_functions import *

#itomic specific
def get_itomic_Data (gene, hub, dataset, samples):
    values= Probe_values(hub, dataset, samples, gene)
    itomic_Data =dict(zip(samples, values))
    return itomic_Data


def Nof1_output (Nof1_sample, original_label, gene, itomic_Data, Nof1_theta):
    print
    print 'Sample: ', Nof1_sample
    print 'Gene:', original_label
    if original_label!= gene:
        print 'HUGO gene name:', gene

    Nof1_value = itomic_Data[Nof1_sample]
    Nof1_TPM = revert_Log2_theta(Nof1_value, Nof1_theta)

    print "log2(TPM):", Nof1_value, "TPM:", '{:.2f}'.format(Nof1_TPM)

def filer_header (comparison_list, fout):
    headerList =["label","gene"]
    header2ndList =["",""]
    for item in comparison_list:
        headerList.append(item["label"]+ ' (n=' + str(len(item["samples"]))+")")
        headerList.append("Range of previous ITOMIC samples vs. " + item["label"])
        header2ndList.append("Rank %")
        header2ndList.append("Rank %")
    headerList.extend(["",""])
    header2ndList.extend(["log2(TPM)","TPM"])

    fout.write(string.join(headerList,'\t') +'\n')
    fout.write(string.join(header2ndList,'\t') +'\n')


def itomic_Nof1(Nof1_item, original_labels, geneMappping, comparison_list, outputfile):
    itomic_samples = dataset_samples(Nof1_item["hub"], Nof1_item["dataset"])

    #file header output
    fout = open(outputfile,'w')
    filer_header (comparison_list, fout)

    for original_label in original_labels:
        gene = geneMappping[original_label];
        itomic_Data = get_itomic_Data (gene, Nof1_item["hub"], Nof1_item["dataset"], itomic_samples)

        #screen output
        Nof1_output (Nof1_item["samples"][0], original_label, gene, itomic_Data, Nof1_item["log2Theta"])


        Nof1_value = itomic_Data[Nof1_item["samples"][0]]

        outputList =[original_label, gene]

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
            outputList.append(string.join(map(lambda x: '100' if (x[1] > 99.5) else '{:.2g}'.format(x[1]),
                r_and_p_values),', '))

            print
            print name +" ( n=", len(h_l_values), "):"
            print "rank:", rank
            print map(lambda x: x[0], r_and_p_values)

            print "expression level percentile:", '{:.2f}%'.format(percentage)

            for list in zip(itomic_Data.keys(), r_and_p_values):
                print list[0], list[1][0], '{:.2f}%'.format(list[1][1])

        outputList.append('{:.2f}'.format(Nof1_value))
        outputList.append('{:.2f}'.format(revert_Log2_theta(Nof1_value, Nof1_item["log2Theta"])))
        fout.write(string.join(outputList,'\t') +'\n')
    fout.close()

def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "expression level percentile is the percentile of samples has lower expression than sample of interest."
    print "High expression level percentile means higher expression."
    print
