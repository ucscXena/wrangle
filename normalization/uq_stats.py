import string, sys, os
import json
import math, numpy

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/../xena/")
import xenaAPI

def getIDs(geneListfile):
    fin = open(geneListfile, 'r')
    IDs = {}
    for line in fin.readlines():
        gene = string.split(string.strip(line),'\t')[0]
        if gene[0] =="#":
            continue
        if gene not in IDs:
            IDs[gene]=0
    fin.close()
    return IDs.keys()

#revert log, include zero, upper quantile + spread
def stats_uq_scale_include_zero_revertLog2 (values, pseudo = 0):
    values = map(lambda x: float(x), values)
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x) else x, values)
    newValues = sorted(values, key=lambda f: float('-inf') if math.isnan(f) else f)
    L = len(values)
    pos = int(L * 0.75)
    uq = newValues[pos]

    # might be single cell data, or data with HUGE drop off, sample should not be normalized
    if uq == 0:
        return None

    values = map(lambda x: math.log((x/uq + pseudo), 2) if not math.isnan(x) else x, values)
    values = filter (lambda x: not math.isnan(x), values)
    var = numpy.var(values)
    sd = math.sqrt(var)
    values.sort()
    L = len(values)

    # might be single cell data, or data with HUGE drop off, sample should not be normalized
    if sd == 0 or (values[int(L * 0.75)]- values[int(L * 0.5)]) == 0:
        return None

    return {
        "uq": uq,
        "log2_uq": math.log(uq,2),
        "log2_sd": sd,
        "log2_75_50": values[int(L * 0.75)]- values[int(L * 0.5)]
    }

def getStats (hub, dataset, samples, mode, pseudo, genes, outputOffset):
    gN = 500
    sN = 100

    #compute uq offset
    params = {}
    for k in range (0, len(samples), sN):
        sList = samples[k:k+sN]
        sample_values =[]
        for i in range (0, len(genes), gN):
            pList = genes[i:i+gN]
            if mode == "probe":
                values = xenaAPI.Probes_values (hub, dataset, sList, pList)
            elif mode == "gene":
                values = xenaAPI.Genes_values (hub, dataset, sList, pList)
            for m in range (0, len(values)):
                if len(values[m]) == 0:
                    values[m] = ['nan'] * len(sList)
            sample_values.extend(values)
            print i
            #if i>gN:
            #    break
        sample_values = zip(*sample_values)

        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            ret = stats_uq_scale_include_zero_revertLog2(values, pseudo)
            if ret:
                params[sample] = ret
                print sample, ret["uq"], ret["log2_uq"], ret["log2_sd"], ret["log2_75_50"]
            else:
                print sample, "with HUGO drop off, should not be normalized"

    #output
    fout_Offset = open(outputOffset, 'w')
    header = params[params.keys()[0]].keys()
    fout_Offset.write("sample\t"+ string.join(header, '\t')+ '\n')
    for sample in samples:
        if sample not in params:
            continue
        sample_param =  params[sample]
        fout_Offset.write(sample)
        for key in header:
            fout_Offset.write('\t'+str(sample_param[key]))
        fout_Offset.write('\n')
    fout_Offset.close()
    return params

if __name__ == "__main__":
    if len(sys.argv[:])!= 4:
        print "python uq_stats.py dataset_json \
            normalization_gene_file(first_column_gene) \
            output_stats\n"
        sys.exit()

    jsonFile = sys.argv[1]
    J = json.loads(open(jsonFile).read())

    hub =J["hub"]
    dataset = J["dataset"]
    if J.has_key("samples"):
        samples = J["samples"]
    else:
        samples = xenaAPI.dataset_samples(hub, dataset)
    mode = J["mode"]
    if J.has_key("pseudo"):
        pseudo = J["pseudo"]
    else:
        pseudo = 0

    genes = getIDs(sys.argv[2])

    outputOffset = sys.argv[3]
    params = getStats(hub, dataset, samples, mode, pseudo, genes, outputOffset)

