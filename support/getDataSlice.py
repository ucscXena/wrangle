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

def process (hub, dataset, samples, mode, genes, outputMatrix_T):
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
                    values[m] = [''] * len(sList)
                    #print pList[m], m, values[m]
            sample_values.extend(values)
            print i
            #if i>gN:
            #    break
        sample_values = zip(*sample_values)
        for j in range(0, len(sList)):
            sample = sList[j]
            values = sample_values[j]
            fout_T.write(sample+'\t')
            fout_T.write(string.join(map(lambda x: str(x), values),'\t')+'\n')

    fout_T.close()


if __name__ == "__main__":
    if len(sys.argv[:])!= 4:
        print "python getDataSlice.py  dataset_json \
            gene_file(\"ALL\" or gene file) \
            outputMatrix_T\n"
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
    if J.has_key("pseudocount"):
        pseudo = J["pseudocount"]
    else:
        print "no pseudocount specified"
        sys.exit()

    genefile = sys.argv[2]
    if genefile == "ALL":
        genes =[]
        mode = "probe"
    else:
        genes = getIDs(genefile)


    outputMatrix_T = sys.argv[3]
    process (hub, dataset, samples, mode, genes, outputMatrix_T)
