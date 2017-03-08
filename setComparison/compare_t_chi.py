import string, sys, os
import json
import scipy.stats
import numpy
import math

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/../xenaAPI/")
import Nof1_functions

def getIDs(geneListfile):
    fin = open(genelistfile, 'r')
    genes = {}
    for line in fin.readlines():
        gene = string.split(string.strip(line),'\t')[0]
        if gene not in genes:
            genes[gene]=0
    fin.close()
    return genes.keys()

#chiq_contigency (add mutual information if the table is 2 x2)
def chi2_contingency_by_probes (hub, dataset, probes, samples1, samples2, output, kp = None):
    fout = open(output,'w')
    o_list =["probe", "chi2_stat", "Gtest_p", "Mutual_information"]
    if kp:
        o_list.append("observed")
        o_list.append("expected")

    fout.write(string.join(o_list, '\t') +'\n')

    N = 100
    for i in range (0, len(probes), N):
        pList = probes[i:i+N]
        values1 = xenaAPI.Probes_values (hub, dataset, samples1, pList)
        values2 = xenaAPI.Probes_values (hub, dataset, samples2, pList)

        for j in range(0, len(pList)):
            probe = pList[j]
            v1 = values1[j]
            v1 = [x for x in v1 if x != 'NaN']

            v2 = values2[j]
            v2 = [x for x in v2 if x != 'NaN']

            codes = numpy.unique(v1 +v2)
            if len(codes) == 1:
                continue

            array1 = map(lambda code : len([x for x in v1 if x == code]), codes)
            array2 = map(lambda code : len([x for x in v2 if x == code]), codes)
            observed = [[array1, array2]]
            try:
                chi2, p, dof, expected = scipy.stats.chi2_contingency(observed, lambda_="log-likelihood")
                o_list = []
                if len(array1) == 2 and len(array2) ==2: #only makes sense for 2 by 2 table, key parameter: 0,1
                    MI = chi2 / (2*len(v1) + len(v2))
                    o_list = [probe, str(chi2), str(p), str(MI)]
                else:
                    o_list = [probe, str(chi2), str(p), '']
                if kp:
                    o_list.append(str(observed[0][kp[0]][kp[1]]))
                    o_list.append(str(expected[0][kp[0]][kp[1]]))
                fout.write(string.join(o_list,'\t')+'\n')
                #print probe, codes, chi2, p, MI, observed, expected
            except:
                print probe, "bad"

    fout.close()

#ttest
def ttest_by_ids (hub1, hub2, dataset1, dataset2, ids, samples1, samples2, mode1, mode2, output):
    fout = open(output,'w')
    o_list =["probe", "t_stat", "test_p", "mean1", "mean2"]
    fout.write(string.join(o_list, '\t') +'\n')

    N = 100
    for i in range (0, len(ids), N):
        pList = ids[i:i+N]
        if mode1 == "probe":
            values1 = xenaAPI.Probes_values (hub1, dataset1, samples1, pList)
        elif mode1 == "gene":
            values1 = xenaAPI.Genes_values (hub1, dataset1, samples1, pList)
        if mode2 == "probe":
            values2 = xenaAPI.Probes_values (hub2, dataset2, samples2, pList)
        elif mode2 == "gene":
            values2 = xenaAPI.Probes_values (hub2, dataset2, samples2, pList)

        for j in range(0, len(pList)):
            probe = pList[j]
            v1 = values1[j]
            v1 = map(lambda x: float(x), v1)
            v1 = [x for x in v1 if math.isnan(x) == 0]
            mean1 = numpy.average (v1)

            v2 = values2[j]
            v2 = map(lambda x: float(x), v2)
            v2 = [x for x in v2 if math.isnan(x) == 0]
            mean2 = numpy.average (v2)

            try:
                tStat, p = scipy.stats.ttest_ind(v1, v2, equal_var=False)
                fout.write(string.join([probe, str(tStat), str(p), str(mean1), str(mean2)], '\t') + '\n')
            except:
                print probe, bad
    fout.close()

if __name__ == "__main__":
    if len(sys.argv[:])!= 5:
        print "python compare_t_chi.py json1 json2 ids(genes/probes) output\n"
        sys.exit()

    jsonFile1 = sys.argv[1]
    J1 = json.loads(open(jsonFile1).read())

    jsonFile2 = sys.argv[2]
    J2 = json.loads(open(jsonFile2).read())

    hub1 =J1["hub"]
    dataset1 = J1["dataset"]
    samples1 = J["samples"]
    mode1 = J1["mode"]

    hub2 =J2["hub"]
    dataset2 = J2["dataset"]
    samples2 = J2["samples"]
    mode2 = J2["mode"]

    ids = getIDs(sys.argv[3])

    output = sys.argv[4]
    ttest_by_probes (hub, dataset_cnv, ids, samples1, samples2, mode1, mode2, output)

    """
    dataset_mut = "mc3.v0.2.8.PUBLIC.nonsilentGene.xena"
    mut_genes = Nof1_functions.dataset_fields(hub, dataset_mut)
    output = 'chisqure_mut'
    key_parameter = [0,1]  #mutants, value =1
    chi2_contingency_by_probes (hub, dataset_mut, mut_genes, samples_mut, samples_wt, output, key_parameter)
    """


