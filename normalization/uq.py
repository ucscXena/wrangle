import string, sys, os
import json
import math, numpy

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/../xena/")
import xenaAPI

uq_target = 32.0
uq_sd_scale_target = 2.5
uq_75_50_scale_target = 1.25

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
def uq_include_zero_revertLog2 (values, parameter, pseudo = 0):
    uq = parameter["uq"]
    #float
    values = map(lambda x: float(x), values)
    #revert log
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x)
        else x, values)
    #upper q
    values = map(lambda x: (math.log((x/uq*uq_target + pseudo), 2)
                if abs(x-pseudo) > 0.001 else math.log(pseudo,2))
            if not math.isnan(x)
            else x, values)
    return values

#revert log, include zero, upper quantile, scale with SD
def uq_SDscale_include_zero_revertLog2 (values, parameter, pseudo = 0):
    uq = parameter["uq"]
    scale = parameter["log2_sd"]
    #float
    values = map(lambda x: float(x), values)
    #revert log
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x)
        else x, values)
    #upper q
    values = map(lambda x: (
            (math.log((x/uq*uq_target + pseudo), 2)- math.log(uq_target + pseudo,2))/scale * uq_sd_scale_target + math.log(uq_target + pseudo,2)
                if abs(x-pseudo) > 0.001 else math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    values = map(lambda x: x if x > math.log(pseudo,2) else (math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    return values

#revert log, include zero, upper quantile, scale with 75-50
def uq_7550scale_include_zero_revertLog2 (values, parameter, pseudo = 0):
    uq = parameter["uq"]
    scale = parameter["log2_75_50"]
    #float
    values = map(lambda x: float(x), values)
    #revert log
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x)
        else x, values)
    #upper q
    values = map(lambda x: (
            (math.log((x/uq*uq_target + pseudo), 2)- math.log(uq_target + pseudo,2))/scale * uq_75_50_scale_target + math.log(uq_target + pseudo,2)
                if abs(x-pseudo) > 0.001 else math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    values = map(lambda x: x if x > math.log(pseudo,2) else (math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    return values

#revert log, include zero, upper quantile + spread
def stats_uq_scale_include_zero_revertLog2 (values, pseudo = 0):
    values = map(lambda x: float(x), values)
    values = map(lambda x: math.pow(2,x) - pseudo if not math.isnan(x) else x, values)
    newValues = sorted(values, key=lambda f: float('-inf') if math.isnan(f) else f)
    L = len(values)
    pos = int(L * 0.75)
    uq = newValues[pos]

    values = map(lambda x: math.log((x/uq*uq_target + pseudo), 2) if not math.isnan(x) else x, values)
    values = filter (lambda x: not math.isnan(x), values)
    var = numpy.var(values)
    sd = math.sqrt(var)
    values.sort()
    L = len(values)
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
            sample_values.extend(values)
            print i
            #if i>500:
            #    break
        sample_values = zip(*sample_values)

        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            ret = stats_uq_scale_include_zero_revertLog2(values, pseudo)
            params[sample] = ret

            print sample, ret["uq"], ret["log2_uq"], ret["log2_sd"], ret["log2_75_50"]

    #output
    fout_Offset = open(outputOffset, 'w')
    header = params[samples[0]].keys()
    fout_Offset.write("sample\t"+ string.join(header, '\t')+ '\n')
    for sample in samples:
        sample_param =  params[sample]
        fout_Offset.write(sample)
        for key in header:
            fout_Offset.write('\t'+str(sample_param[key]))
        fout_Offset.write('\n')
    fout_Offset.close()
    return params

def process (hub, dataset, samples, mode, pseudo, method, outputMatrix_T, params):
    fout_T = open(outputMatrix_T, 'w')

    gN = 500
    sN = 100

    #convert data
    probes = xenaAPI.dataset_fields(hub, dataset)
    probes.remove("sampleID")
    fout_T.write('sample\t'+string.join(probes,"\t")+'\n')

    for k in range (0, len(samples), sN):
        sList = samples[k:k+sN]
        sample_values =[]
        for i in range (0, len(probes), gN):
            pList = probes[i:i+gN]
            values = xenaAPI.Probes_values (hub, dataset, sList, pList)
            sample_values.extend(values)
            print i
            #if i>500:
            #    break
        sample_values = zip(*sample_values)

        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            parameter = params[sample]

            values = method(values, parameter, pseudo)

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

    params = getStats(hub, dataset, samples, mode, pseudo, genes, outputOffset)

    #method = uq_include_zero_revertLog2
    #method = uq_SDscale_include_zero_revertLog2
    method = uq_7550scale_include_zero_revertLog2

    process (hub, dataset, samples, mode, pseudo, method, outputMatrix_T, params)
