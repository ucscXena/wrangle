import string, sys, os
import json
import math, numpy

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/../xena/")
import xenaAPI

uq_target = 28.87894851
uq_scale_target = 2.336337271

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

#revert log, include zero, upper quantile
def uq_include_zero_revertLog2 (values, pseudo = 0):
    values = map(lambda x: float(x), values)
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x) else x, values)
    newValues = sorted(values, key=lambda f: float('-inf') if math.isnan(f) else f)
    L = len(values)
    pos = int(L * 0.75)
    uq = newValues[pos]
    return {
        "uq": uq,
        "log2_uq": math.log(uq,2),
        "log2_scale" : uq_scale_target
    }

#revert log, include zero, upper quantile + spread
def uq_scale_include_zero_revertLog2 (values, pseudo = 0):
    values = map(lambda x: float(x), values)
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x) else x, values)
    newValues = sorted(values, key=lambda f: float('-inf') if math.isnan(f) else f)
    L = len(values)
    pos = int(L * 0.75)
    uq = newValues[pos]

    values = map(lambda x: math.log((x/uq + pseudo), 2) if not math.isnan(x) else x, values)
    values = filter (lambda x: not math.isnan(x), values)
    var = numpy.var(values)
    sd = math.sqrt(var)
    values.sort()
    L = len(values)
    return {
        "uq": uq,
        "log2_uq": math.log(uq,2),
        "log2_scale": sd,
        "log2_75_50": values[int(L * 0.75)]- values[int(L * 0.5)]
    }

def process (hub, dataset, samples, mode, pseudo, genes, method,
    outputMatrix_T, outputOffset):

    fout_T = open(outputMatrix_T, 'w')
    fout_Offset = open(outputOffset, 'w')

    gN = 500
    sN = 100

    fout_Offset.write(string.join(["sample", "uq", "log2_uq", "log2_scale", "log2_75_50"], '\t')+ '\n')

    #compute uq offset
    offsets = {}
    log2scales = {}
    for k in range (0, len(samples), sN):
        sList = samples[k:k+sN]
        sample_values =[]
        for i in range (0, len(genes), gN):
            pList = genes[i:i+gN]
            if mode == "probe":
                values = xenaAPI.Probes_values (hub, dataset, sList, pList)
            elif mode == "gene":
                values = xenaAPI.Genes_values (hub, dataset, sList, pList)
            sample_values.extend(values)
            print i

        sample_values = zip(*sample_values)

        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            ret = method(values, pseudo)
            offsets[sample] = ret["uq"]
            log2scales[sample] = ret["log2_scale"]

            print sample, ret["uq"], ret["log2_uq"], ret["log2_scale"]
            #output
            fout_Offset.write(
                string.join([sample, str(ret["uq"]), str(ret["log2_uq"]), str(ret["log2_scale"]), str(ret["log2_75_50"])], '\t')
                + '\n')

    fout_Offset.close()

    #convert data
    probes = xenaAPI.dataset_fields(hub, dataset)
    fout_T.write('sample\t'+string.join(probes,"\t")+'\n')

    for k in range (0, len(samples), sN):
        sList = samples[k:k+sN]
        sample_values =[]
        for i in range (0, len(probes), gN):
            pList = probes[i:i+gN]
            values = xenaAPI.Probes_values (hub, dataset, sList, pList)
            sample_values.extend(values)
            print i

        sample_values = zip(*sample_values)

        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            uq = offsets[sample]
            log2_scale = log2scales[sample]
            values = map(lambda x: float(x), values)
            values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x) else x, values)
            values = map(lambda x: (math.log((x/uq*uq_target + pseudo), 2)- math.log(uq_target + pseudo,2))/log2_scale * uq_scale_target if not math.isnan(x) else x, values)

            fout_T.write(sample+'\t')
            fout_T.write(string.join(map(lambda x: str(x), values),'\t')+'\n')

    fout_T.close()


if __name__ == "__main__":
    if len(sys.argv[:])!= 5:
        print "python uq.py dataset_json gene_file(first_column_gene) outputMatrix_T outputOffset\n"
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

    outputMatrix_T = sys.argv[3]
    outputOffset = sys.argv[4]

    process (hub, dataset, samples, mode, pseudo, genes, uq_scale_include_zero_revertLog2,
        outputMatrix_T, outputOffset)
