import string, sys
import json
import xena_query as xena

if len(sys.argv[:]) < 4:
    print "python validateDefaultDataset.py input_default_txt type(txt/json) huburl(s)"
    sys.exit()

hubs = sys.argv[3:]
type = sys.argv[2]

if type not in ["txt","json"]:
    print "bad type"
    sys.exit()

input = sys.argv[1]
fin = open(input,'U')

if type == "txt":
    fin.readline()
    for line in fin.readlines():
        list = string.split(string.strip(line),'\t')
        #print list
        cohort = list[0]
        datasets = list[1:]
        db_d_list = map(lambda hub: xena.datasets_list_in_cohort(hub, cohort), hubs)
        reduce(lambda x, y: x.extend(y), db_d_list)
        db_d_list = db_d_list[0]
        for d in datasets:
            if d =="":
                continue
            if d not in db_d_list:
                print "not found", cohort,d

elif type == "json":
    J = json.load(fin)
    for cohort in J.keys():
        db_d_list = map(lambda hub: xena.datasets_list_in_cohort(hub, cohort), hubs)
        reduce(lambda x, y: x.extend(y), db_d_list)
        db_d_list = db_d_list[0]

        for dataset in J[cohort].values():
            if dataset == "":
                continue
            if dataset not in db_d_list:
                print "not found", cohort, dataset

fin.close()
