import string, math
import xena_query as xena

def rank_and_percentage(v, v_list):
    v = float(v)
    l = map(float, v_list)
    l.sort()
    l.reverse() #large to small
    K = len(l)
    for i in range(0, K):
        if v > l[i]:
            return i, float(i)/K *100
    return K, 100.0

def Gene_values(hub, dataset, samples, gene):
    values = xena.dataset_gene_values (hub, dataset, samples, [gene])
    return values[0]["scores"][0]

def Probe_values(hub, dataset, samples, probe):
    values = xena.dataset_probe_values (hub, dataset, samples, [probe])
    return values[0]

def dataset_samples(hub,dataset):
    return xena.dataset_samples(hub, dataset)



