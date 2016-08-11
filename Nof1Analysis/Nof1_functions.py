import string, math
import xena_query as xena

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
            return i+1, float(i)/K *100
    return K+1, 100.0

def checkSamples (Nof1_sample, hub, dataset):
    samples = dataset_samples(hub, dataset)

    if Nof1_sample not in samples:
        print Nof1_sample, "is not in dataset \""+ dataset +"\"."
        print
        print "Find out the dataset samples:", "https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host="+ hub + "&dataset="+ dataset + "&allSamples=true"
        print
        return 1
    else:
        return 0

def checkFields (fields, mapping, hub, dataset):
    d_fields = dataset_fields (hub, dataset)
    for field in fields:
        if field  in d_fields:
           mapping[field] = field
        elif string.upper(field) in d_fields:
            mapping[field] = string.upper(field)
        elif field in mapping and mapping[field] in d_fields:
            pass
        else:
            print field, "is not in dataset \""+ dataset +"\"."
            print
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



