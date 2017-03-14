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

def getParams(paramfile):
    fin = open(paramfile, 'r')
    params = {}
    data = string.split(string.strip(fin.readline()),'\t')
    header = {}

    for i in range(1, len(data)):
        header[i] = data[i]
    for line in fin.readlines():
        data = string.split(string.strip(line),'\t')
        sample = data[0]
        if sample not in params:
            params[sample]={}
        for i in range(1, len(data)):
            key = header[i]
            params[sample][key] = float(data[i])

    fin.close()
    return params

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
                if x-pseudo > 0.001 else math.log(pseudo,2))
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
                if x-pseudo > 0.001 else math.log(pseudo,2)) if not math.isnan(x)
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
                if x-pseudo > 0.001 else math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    values = map(lambda x: x if x > math.log(pseudo,2) else (math.log(pseudo,2)) if not math.isnan(x)
            else x, values)
    return values

def process (hub, dataset, samples, mode, pseudo, method, genes, outputMatrix_T, params):
    fout_T = open(outputMatrix_T, 'w')

    gN = 500
    sN = 100

    #convert data
    if len(genes)==0:
        probes = xenaAPI.dataset_fields(hub, dataset)
        probes.remove("sampleID")
    else:
        probes = genes
    fout_T.write('sample\t'+string.join(probes,"\t")+'\n')

    for k in range (0, len(samples), sN):
        sList = samples[k:k+sN]
        sample_values =[]
        for i in range (0, len(probes), gN):
            pList = probes[i:i+gN]
            if mode =="probe":
                values = xenaAPI.Probes_values (hub, dataset, sList, pList)
            else:
                values = xenaAPI.Genes_values (hub, dataset, sList, pList)
            for m in range (0, len(values)):
                if len(values[m]) == 0:
                    values[m] = ['nan'] * len(sList)
                    #print pList[m], m, values[m]
            sample_values.extend(values)
            print i
            #if i>gN:
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
    if len(sys.argv[:])!= 6:
        print "python uq.py  dataset_json \
            whole_genome_gene_file(\"ALL\" or gene file) \
            parameter_file \
            method \
            outputMatrix_T\n"
        print "method: uq_include_zero_revertLog2"
        print "        uq_SDscale_include_zero_revertLog2"
        print "        uq_7550scale_include_zero_revertLog2"
        print
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

    genefile = sys.argv[2]
    if genefile == "ALL":
        genes =[]
        mode = "probe"
    else:
        genes = getIDs(genefile)

    #paramters read from file
    params = getParams(sys.argv[3])

    methodname = sys.argv[4]
    if methodname not in [
        "uq_include_zero_revertLog2",
        "uq_SDscale_include_zero_revertLog2",
        "uq_7550scale_include_zero_revertLog2"
        ]:
        print "wrong method name\n"
        sys.exit()
    else:
        if methodname == "uq_include_zero_revertLog2":
            method = uq_include_zero_revertLog2
        elif methodname ==  "uq_SDscale_include_zero_revertLog2":
            method =  uq_SDscale_include_zero_revertLog2
        elif methodname ==    "uq_7550scale_include_zero_revertLog2":
            method = uq_7550scale_include_zero_revertLog2

    outputMatrix_T = sys.argv[5]
    process (hub, dataset, samples, mode, pseudo, method, genes, outputMatrix_T, params)
