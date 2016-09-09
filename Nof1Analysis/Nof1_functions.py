import string, math
import xena_query as xena

def average (v_list):
    return reduce(lambda x, y: x+y, v_list)/float(len(v_list))

def standard_deviation (v_list):
    m = average (v_list)
    variance = reduce(lambda x, y : x + (y -m) *(y-m), v_list)/float(len(v_list))
    SD = math.sqrt(variance)
    return SD

def clean_then_sort_high_to_low (v_list):
    l = filter (lambda x: not math.isnan(x), map(float, v_list))
    l.sort()
    l.reverse() #large to small
    return l

def rank_and_percentage(v, v_list):
    v = float(v)
    K = len(v_list)
    for i in range(0, K):
        if v > v_list[i]:
            return i+1, 100 - float(i)/K *100
    return K+1, 0.0

def revert_Log2_theta (Nof1_value, theta):
    return math.pow(2, Nof1_value)- theta

def checkSamples (Nof1_sample, hub, dataset):
    samples = dataset_samples(hub, dataset)

    if Nof1_sample not in samples:
        print Nof1_sample, "is not found in dataset."
        print
        print "Find out the dataset samples:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host="+ hub + "&dataset="+ dataset + "&allSamples=true&label="+dataset
        print
        return 1
    else:
        return 0

def checkFields (fields, mapping, hub, dataset, cleanFuntion):
    d_fields = dataset_fields (hub, dataset)
    for original_field in fields:
        field = cleanFuntion (original_field)
        if field in d_fields:
           mapping[original_field] = field
        elif string.upper(field) in d_fields:
            mapping[original_field] = string.upper(field)
        elif original_field in mapping and mapping[original_field] in d_fields:
            pass
        else:
            print original_field, "is not found in dataset."
            print
            print "Find out the dataset identifiers:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host="+ hub + "&dataset="+ dataset + "&allFields=true"
            return 1
    return 0

def Gene_values (hub, dataset, samples, gene):
    values = xena.dataset_gene_values (hub, dataset, samples, [gene])
    return values[0]["scores"][0]

def Probe_values (hub, dataset, samples, probe):
    values = xena.dataset_probe_values (hub, dataset, samples, [probe])
    return values[0]

def dataset_samples (hub,dataset):
    return xena.dataset_samples(hub, dataset)

def dataset_fields (hub, dataset):
    return xena.dataset_field (hub, dataset)



