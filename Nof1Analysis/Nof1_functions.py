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

def Gene_values(hub, dataset, samples, gene):
    values = xena.dataset_gene_values (hub, dataset, samples, [gene])
    return values[0]["scores"][0]

def Probe_values(hub, dataset, samples, probe):
    values = xena.dataset_probe_values (hub, dataset, samples, [probe])
    return values[0]

def dataset_samples(hub,dataset):
    return xena.dataset_samples(hub, dataset)



