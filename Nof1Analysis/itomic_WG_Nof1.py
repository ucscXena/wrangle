import string, sys
import scipy.stats
import numpy
from Nof1_functions import *

#itomic specific
def get_itomic_Data (genes, hub, dataset, samples):
    dic ={}
    values_list= Probes_values(hub, dataset, samples, genes)
    ret_list = map(lambda values: dict(zip(samples, values)), values_list)
    return ret_list  # a list of dictionaries, orders by genes

def getMatrixData(file):
    fin = open(file,'r')
    line = fin.readline()
    dic={}
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line[:-1],'\t')
        id = data [0]
        dic [id] = data[1:]
    fin.close()
    return dic

def file_header (comparison_item, Nof1_item, fout):
    topline = "ITOMIC samples vs. "+ comparison_item["label"] + ' (n=' + str(len(comparison_item["samples"])) + ")"
    headerList =["label","gene"]
    headerList.append("logTPM p")
    headerList.append("logTPM t")
    headerList.append("logTPM ave itomic")
    headerList.append("logTPM ave2")
    headerList.append("rank % p")
    headerList.append("rank % t")
    headerList.append("rank % ave")
    headerList.append("rank % SD")

    for sample in Nof1_item["samples"]:
        headerList.append("Rank %")

    fout.write(topline+'\n')
    fout.write(string.join(headerList,'\t') +'\n')


def itomic_Nof1(Nof1_item, original_labels, geneMappping, comparison_item, outputfile):
    itomic_samples = dataset_samples(Nof1_item["hub"], Nof1_item["dataset"])

    fout = open(outputfile,'w')
    foutdata = open(outputfile+"_data",'w') #pure data file

    #full file header output
    file_header (comparison_item, Nof1_item, fout)

    #data file header
    foutdata.write("gene")
    foutdata.write('\t'+ string.join(Nof1_item["samples"],'\t'))
    foutdata.write('\t'+ string.join(Nof1_item["samples"],'\t'))
    foutdata.write('\t'+ string.join(Nof1_item["samples"],'\t'))
    foutdata.write('\n')

    foutdata.write("gene")
    foutdata.write('\t'+ string.join(map(lambda x : "RANK %", Nof1_item["samples"]),'\t'))
    foutdata.write('\t'+ string.join(map(lambda x : "Log2TPM", Nof1_item["samples"]),'\t'))
    foutdata.write('\t'+ string.join(map(lambda x : "TPM", Nof1_item["samples"]),'\t'))
    foutdata.write('\n')


    # comparison data
    if "file" in comparison_item:
        file = comparison_item["file"]
        cData = getMatrixData(file)
    else:
        cData = None

    hub = comparison_item["hub"]
    dataset = comparison_item["dataset"]
    samples = comparison_item["samples"]
    name = comparison_item["name"]
    mode = comparison_item["mode"]

    if cData:
        n = 2000
    else:
        n = int(100000/len(samples))
        if n < 100:
            n =100
        print n

    for k in range (0, len(original_labels), n):
        labels = original_labels[k:k+n]
        genes = map(lambda original_label: geneMappping[original_label] if (original_label in geneMappping) else original_label,labels)

        print genes[:10], "..."

        #all itomic data
        all_data_list = get_itomic_Data (genes, Nof1_item["hub"], Nof1_item["dataset"], itomic_samples)

        if cData:
            compare_data_list = None
        else:
            # get data for comparison
            if mode == "gene":
                compare_data_list = Genes_values (hub, dataset, samples, genes)
            if mode == "probe":
                compare_data_list = Probes_values (hub, dataset, samples, genes)

        for m in range (0, n):
            if m  == len(genes):
                break
            gene = genes[m]
            label = labels[m]
            outputList =[label, gene]
            data_outputList = [gene]
            all_Data = all_data_list[m]

            itomic_values = map(lambda sample: all_Data[sample], Nof1_item["samples"])

            #all itomic
            allsample_Data = all_data_list[m]

            # cohort data
            if cData:
                if gene in cData:
                    values = cData[gene]
                else:
                    values =[]
                #print len(values)
                #h_l_values = clean (values)
                #r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), itomic_values)
                #print r_and_p_values

            if compare_data_list:
                compare_gene_obj = compare_data_list[m]
                values = compare_gene_obj['scores'][0] #############
                #print len(values)
                #h_l_values = clean (values)
                #r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), itomic_values)
                #print r_and_p_values

            h_l_values = clean (values)
            r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), itomic_values)

            #ttest p value
            try:
                tStat, p = scipy.stats.ttest_ind(allsample_Data.values(), h_l_values, equal_var=False)
                mean1 = numpy.mean( allsample_Data.values())
                mean2 = numpy.mean( h_l_values)
                outputList.append ('{:.4f}'.format(p)) # ttest p value
                outputList.append (str(tStat)) # ttest t
                outputList.append (str(mean1))
                outputList.append (str(mean2))
            except:
                outputList.append ('')
                outputList.append ('')
                outputList.append ('')

            #SD
            all_r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), allsample_Data.values())
            SD = standard_deviation(map(lambda x: x[1], all_r_and_p_values))

            try:
                r_list = map(lambda x: x[1], r_and_p_values)
                tStat, p = scipy.stats.ttest_ind(r_list, range(0,100), equal_var=False)
                mean1 = numpy.mean(r_list)
                outputList.append ('{:.4f}'.format(p)) # ttest p value
                outputList.append (str(tStat)) # ttest t
                outputList.append (str(mean1))
            except:
                outputList.append ('')
                outputList.append ('')
                outputList.append ('')

            outputList.append ('{:.2f}'.format(SD))#rank SD
            outputList.extend(map(lambda x: '{:.2f}%'.format(x[1]), r_and_p_values)) #rank %

            data_outputList.extend(map(lambda x: '{:.2f}%'.format(x[1]), r_and_p_values)) #rank %
            data_outputList.extend(map(lambda x: '{:.2f}'.format(x), itomic_values)) #Log2TPM
            data_outputList.extend(map(lambda x: '{:.2f}'.format(revert_Log2_theta(x, Nof1_item["log2Theta"])), itomic_values)) #TPM

            fout.write(string.join(outputList,'\t') +'\n')
            foutdata.write(string.join(data_outputList,'\t') +'\n')

    fout.write("\n")
    fout.write("Rank % : percentile of samples with lower expression than sample of interest.\n")
    fout.write("Higher Rank %  means higher expression.\n")
    fout.close()
    foutdata.close()

def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "Rank % is the percentile of samples with lower expression than sample of interest."
    print "Higher Rank %  means higher expression."
    print
