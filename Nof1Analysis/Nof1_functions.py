import string, math, sys
import xenaPython as xena

def average (v_list):
    return reduce(lambda x, y: x+y, v_list, 0)/float(len(v_list))

def standard_deviation (v_list):
    m = average (v_list)
    variance = reduce(lambda x, y : x + (y -m) *(y-m), v_list, 0)/float(len(v_list))
    SD = math.sqrt(variance)
    return SD

#def clean_then_sort_high_to_low (v_list):
#    l = filter (lambda x: not math.isnan(x), map(float, v_list))
#    l.sort()
#    l.reverse() #large to small
#    return l

def clean (v_list):
    l = filter (lambda x: not math.isnan(x), map(float, v_list))
    return l

#old
#def rank_and_percentage (v, v_list):
#    v = float(v)
#    K = len(v_list)
#    for i in range(0, K):
#        if v > v_list[i]:
#            return i+1, 100 - float(i)/K *100
#    return K+1, 0.0

def rank_and_percentage (v, v_list):
    v_list.sort()  #low to high
    v = float(v)
    K = len(v_list)
    for i in range(0, K):
        if v < v_list[i]:
            return K - i +1, i * 100.0/K
    return 1, 100.0

def revert_Log2_theta (Nof1_value, theta):
    return math.pow(2, Nof1_value)- theta

def checkSamples (Nof1_sample, hub, dataset):
    samples = xena.xenaAPI.dataset_samples(hub, dataset)

    if Nof1_sample not in samples:
        print Nof1_sample, "is not found in dataset."
        print
        print "Find out the dataset samples:", "https://xenabrowser.net/datapages/?host="+ hub + "&dataset="+ dataset + "&allSamples=true&label="+dataset
        print
        return 1
    else:
        return 0

def checkFields (fields, mapping, hub, dataset, cleanFuntion):
    d_fields = xena.xenaAPI.dataset_fields (hub, dataset)
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
            print "Find out the dataset identifiers:", "https://xenabrowser.net/datapages/?host="+ hub + "&dataset="+ dataset + "&allIdentifiers=true"
            return 1
    return 0




