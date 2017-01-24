import string, sys
import json
import xena_query as xena

if len(sys.argv[:]) != 2:
    print "python validateDefaultDataset.py metadata_json"
    sys.exit()


input = sys.argv[1]
fin = open(input,'U')

J = json.load(fin)
for cohort in J.keys():
    db_d_dic ={}
    for dataSubType in J[cohort].keys():
        hub = J[cohort][dataSubType]["host"]
        dataset = J[cohort][dataSubType]["dataset"]
        if hub not in db_d_dic:
            db_d_dic[hub] =  xena.datasets_list_in_cohort(hub, cohort)

        if dataset not in db_d_dic[hub]:
            print "not found", cohort, hub, dataset

fin.close()
