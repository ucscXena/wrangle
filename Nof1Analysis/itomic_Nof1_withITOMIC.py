import string, sys
import scipy.stats
import numpy
import statsmodels.sandbox.stats.multicomp
from Nof1_functions import *

#itomic specific
def get_itomic_Data (gene, hub, dataset, samples):
    values= xenaAPI.Probe_values(hub, dataset, samples, gene)
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

def filer_header (comparison_list, itomic_samples, Nof1_item, fout):
    headerList =["label","gene"]
    header2ndList =["",""]
    headerList.append(Nof1_item["label"]+ ' (n=' + str(len(itomic_samples))+")")
    header2ndList.append("Rank %")

    for item in comparison_list:
        headerList.append(Nof1_item["label"]+ " samples (n=" + str(len(itomic_samples)) + ')' +
            "vs. " +
            item["label"] + ' (n=' + str(len(item["samples"]))+")") #p
        headerList.append("") # t
        headerList.append("") # mean itomic
        headerList.append("") # mean 2

        header2ndList.append("ttest p")
        header2ndList.append("ttest t")
        header2ndList.append("mean itomic")
        header2ndList.append("mean " + item["label"])

    Nof1_sample = Nof1_item["samples"][0]
    headerList.extend([Nof1_sample, Nof1_sample])
    header2ndList.extend(["log2(TPM)","TPM"])

    fout.write(string.join(headerList,'\t') +'\n')
    fout.write(string.join(header2ndList,'\t') +'\n')


def itomic_Nof1(Nof1_item, original_labels, geneMappping, external_comparison_list, outputfile):
    itomic_samples = xenaAPI.dataset_samples(Nof1_item["hub"], Nof1_item["dataset"])
    Nof1_sample = Nof1_item["samples"][0]

    #file header output
    fout = open(outputfile,'w')
    filer_header (external_comparison_list, itomic_samples, Nof1_item, fout)

    pDic = {} # pvalue collection
    for original_label in original_labels:
        if original_label in geneMappping:
            gene = geneMappping[original_label]
        else:
            gene = original_label
        outputList =[original_label, gene]

        itomic_Data = get_itomic_Data (gene, Nof1_item["hub"], Nof1_item["dataset"], itomic_samples)

        #screen output
        Nof1_output (Nof1_sample, original_label, gene, itomic_Data, Nof1_item["log2Theta"])

        Nof1_value = itomic_Data[Nof1_sample]

        rank, percentage =  rank_and_percentage (Nof1_value, itomic_Data.values())
        outputList.append('{:.2f}%'.format(percentage))


        for item in external_comparison_list:
            hub = item["hub"]
            dataset = item["dataset"]
            samples = item["samples"]
            name = item["name"]
            mode = item["mode"]
            gene = string.upper(gene)

            if mode == "gene":
                values = xenaAPI.Gene_values(hub, dataset, samples, gene)
            elif mode == "probe":
                values = xenaAPI.Probe_values(hub, dataset, samples, gene)
            h_l_values = clean (values) # comparison larger cohort values

            # log2TPM statistics
            # ttest p value
            try:
                tStat, p = scipy.stats.ttest_ind(itomic_Data.values(), h_l_values, equal_var=False)
                mean1 = numpy.mean( itomic_Data.values())
                mean2 = numpy.mean( h_l_values)
                outputList.append (str(p)) # ttest p value
                outputList.append (str(tStat)) # ttest t
                outputList.append (str(mean1))
                outputList.append (str(mean2))
                pDic[gene] = p
            except:
                outputList.append('')
                outputList.append('')
                outputList.append('')
                outputList.append('')

            print
            print name +" ( n=", len(h_l_values), "):"
            print "rank:", rank
            #print map(lambda x: x[0], r_and_p_values)

            print "Rank %:", '{:.2f}%'.format(percentage)

            #for list in zip(itomic_Data.keys(), r_and_p_values):
            #    print list[0], list[1][0], '{:.2f}%'.format(list[1][1])

        outputList.append('{:.2f}'.format(Nof1_value))
        outputList.append('{:.2f}'.format(revert_Log2_theta(Nof1_value, Nof1_item["log2Theta"])))
        fout.write(string.join(outputList,'\t') +'\n')
    fout.write("\n")
    fout.write("Rank % : percentile of samples with lower expression than sample of interest.\n")
    fout.write("Higher Rank %  means higher expression.\n")
    fout.close()

def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "Rank % is the percentile of samples with lower expression than sample of interest."
    print "Higher Rank %  means higher expression."
    print
