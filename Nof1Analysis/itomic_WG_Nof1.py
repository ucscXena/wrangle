import string, sys, os
import uuid
import scipy.stats
import numpy
import statsmodels.sandbox.stats.multicomp
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
    headerList.append("logTPM fdr_bh p")
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

    tmpfile = str(uuid.uuid4())
    fout = open(tmpfile,'w')
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


    pDic = {} # all p values for multiple hypo adjustment

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

            #all itomic
            allsample_Data = all_data_list[m]

            # cohort data
            if cData:
                if gene in cData:
                    values = cData[gene]
                else:
                    values =[]

            if compare_data_list:
                compare_gene_obj = compare_data_list[m]
                values = compare_gene_obj['scores'][0] #############

            h_l_values = clean (values)

            if len(h_l_values) == 0: #no comparison data
                continue

            #ttest p value
            try:
                tStat, p = scipy.stats.ttest_ind(allsample_Data.values(), h_l_values, equal_var=False)
                mean1 = numpy.mean( allsample_Data.values())
                mean2 = numpy.mean( h_l_values)
                outputList.append (str(p)) # ttest p value
                outputList.append ('')
                outputList.append (str(tStat)) # ttest t
                outputList.append (str(mean1))
                outputList.append (str(mean2))
                pDic[gene] = p
            except:
                continue

            #SD
            all_r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), allsample_Data.values())
            SD = standard_deviation(map(lambda x: x[1], all_r_and_p_values))

            try:
                r_list = map(lambda x: x[1], all_r_and_p_values)
                tStat, p = scipy.stats.ttest_ind(r_list, range(0,100), equal_var=False)
                mean1 = numpy.mean(r_list)
                outputList.append (str(p)) # ttest p value
                outputList.append (str(tStat)) # ttest t
                outputList.append (str(mean1))
                outputList.append ('{:.2f}'.format(SD))#rank SD
            except:
                continue


            # per sample data output
            itomic_values = map(lambda sample: all_Data[sample], Nof1_item["samples"])
            r_and_p_values = map(lambda x: rank_and_percentage(x, h_l_values), itomic_values)
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

    #add multiple hypo adjusted p values to file
    #http://statsmodels.sourceforge.net/devel/generated/statsmodels.sandbox.stats.multicomp.fdrcorrection0.html
    pCorDic = {}
    genes = pDic.keys()
    rejected, pvalue_corrected =  statsmodels.sandbox.stats.multicomp.fdrcorrection0( map( lambda x : pDic[x], genes),
        alpha=0.05, method='indep', is_sorted=False)
    for i in range(0, len(genes)):
        gene = genes[i]
        pCorDic[gene] = pvalue_corrected[i]

    fout = open(outputfile, 'w')
    fin = open(tmpfile, 'r')
    fout.write(fin.readline())
    fout.write(fin.readline())
    while 1:
        line = fin.readline()
        if line == '' or line == '\n':
            break
        data = string.split(line, '\t')
        gene = data[0]
        data[3] = str(pCorDic[gene])
        fout.write(string.join(data,'\t'))
    fin.close()
    fout.close()
    os.system("rm -f " + tmpfile)

def itomic_legend():
    print "\nExpression values are sorted from high to low."
    print "Low rank means high expression."
    print "Rank % is the percentile of samples with lower expression than sample of interest."
    print "Higher Rank %  means higher expression."
    print
