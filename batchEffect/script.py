import string, math, sys, os
import scipy.stats
import numpy

sys.path.insert(0,os.path.dirname(os.path.realpath(__file__))+"/../xena/")
import xena_datasetlist, gene_list
import xenaAPI

#coding_genes =["APOOP5","ARGFX","ESR1","blah"]
coding_genes = gene_list.protein_coding_genes["genes"]

def clean (v_list):
    l = filter (lambda x: not math.isnan(x), map(float, v_list))
    return l

def parseUnit(unit):  # working progress
    if string.find(unit,"log") != -1:
        logBase = string.strip(string.split(unit,'(')[0])[3:]
        theta = string.strip(string.split(string.split(unit,'(')[1], ')')[0])
        theta = string.strip(string.split(theta, '+')[1])
        if logBase != '' and theta != '':
            return [int(logBase), float(theta)]
    else:
        return None

def mean_variance (values, unit = None):
    cValues = clean(values)
    mean = numpy.mean(cValues)
    var = numpy.var(cValues)
    return {
        "n": len(cValues),
        "mean": mean,
        "var": var
    }

def mean_variance_countZero_countOne (values, unit):
    r = parseUnit(unit)
    if r == None:
        return None
    logBase, theta = r
    zeroInLog = math.log(theta, logBase)
    oneInLog = math.log(1+theta, logBase)

    cValues = clean(values)
    mean = numpy.mean(cValues)
    var = numpy.var(cValues)
    zeros = filter(lambda x : abs(x - zeroInLog) < 0.001, cValues)
    ones = filter(lambda x : (x - oneInLog) < 0.001, cValues)
    return {
        "n": len(cValues),
        "mean": mean,
        "var": var,
        "zeros": len(zeros),
        "ones": len(ones)
    }

def mean_variance_countZero_countOne_countPointOne (values, unit):
    r = parseUnit(unit)
    if r == None:
        return None
    logBase, theta = r
    zeroInLog = math.log(theta, logBase)
    pointOneInLog = math.log(0.1+theta, logBase)
    oneInLog = math.log(1+theta, logBase)

    cValues = clean(values)
    mean = numpy.mean(cValues)
    var = numpy.var(cValues)
    zeros = filter(lambda x : abs(x - zeroInLog) < 0.001, cValues)
    pointOnes = filter(lambda x : (x - pointOneInLog) < 0.001, cValues)
    ones = filter(lambda x : (x - oneInLog) < 0.001, cValues)

    return {
        "n": len(cValues),
        "mean": mean,
        "var": var,
        "zeros": len(zeros),
        "pointOnes": len(pointOnes),
        "ones": len(ones)
    }

def revertLog2_mean_variance_countZero_countOne (values, unit):  ############ need test
    r = parseUnit(unit)
    if r == None:
        return None
    logBase, theta = r
    cValues = clean(values)
    cValues = map(lambda x: math.pow(logBase, x) - theta, cValues)
    mean = numpy.mean(cValues)
    var = numpy.var(cValues)
    zeros = filter(lambda x : abs(x) < 0.001, cValues)
    ones = filter(lambda x : x < 1, cValues)

    return {
        "n": len(cValues),
        "mean": mean,
        "var": var,
        "zeros": len(zeros),
        "ones": len(ones)
    }

def process (obj, IDs, outfile, action):
    fout = open(outfile,'w')
    n = int(100000/len(obj["samples"]))
    if n < 20:
        n = 20
    if n > 500:
        n = 500
    print n

    header = True
    for k in range (0, len(IDs), n):
        ids = IDs[k:k+n]
        if obj['mode'] == "gene":
            values_list = xenaAPI.Genes_values(obj['hub'], obj['dataset'], obj['samples'], ids)
        elif obj['mode'] == "probe":
            values_list = xenaAPI.Probes_values(obj['hub'], obj['dataset'], obj['samples'], ids)

        for i in range(0, len(ids)):
            id = ids[i]
            if obj['mode'] == "gene":
                values = values_list[i]["scores"][0]
            elif obj['mode'] == "probe":
                values = values_list[i]

            ret = action(values, obj["unit"])
            output(id, ret, fout, header)
            if header:
                header = False
    fout.close()


def output(id, stats_obj, fout, header = False):
    if header:
        fout.write('id\t' + string.join(stats_obj.keys(),'\t') + '\n')
    fout.write(id+'\t')
    fout.write(string.join(map(lambda x : str(x), stats_obj.values()),'\t')+'\n')
    print id, stats_obj



if len(sys.argv[:])!= 3:
    print "python script.py dataset_obj_file(json) output_stats"
    sys.exit()

inputfile = sys.argv[1]
fin = open(inputfile,'r')
obj = json.load(fin)
fin.close()

outfile = sys.argv[2]
# obj = xena_datasetlist.TCGA_tumors_geneExp
#outfile = 'TCGA_tumors_log2TPM_stats'

process (obj, coding_genes, outfile, mean_variance_countZero_countOne_countPointOne)

#outfile = 'TCGA_breasts_log2TPM'
#process (xena_datasetlist.TCGA_BRCA_tumors_geneExp, coding_genes, outfile, mean_variance_countZero_countOne)

#outfile = 'itomic_log2TPM'
#process (xena_datasetlist.itomic_geneExp, coding_genes, outfile, mean_variance_countZero_countOne)

