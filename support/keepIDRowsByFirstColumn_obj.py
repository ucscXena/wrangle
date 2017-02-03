import string, sys, os
import json
sys.path.insert(0,"../xena/")
import xenaAPI

def listing (mfile):
    fin = open(mfile,'r')
    dic= {}
    for line in fin.readlines():
        data = string.strip(string.split(line,'\t')[0])
        if data not in dic:
            dic[data]=0
        else:
            print data
    fin.close()
    return dic

def keepProbes (obj, outputFile, keep_dic):
    fout = open(outputFile,'w')
    n = int(100000/len(obj["samples"]))
    if n < 20:
        n = 20
    if n > 500:
        n = 500
    print n

    fout.write("id\t" + string.join(obj['samples'], '\t') + '\n')

    IDs = keep_dic.keys()
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
            fout.write(id + map(lambda x : str(x), values) +'\n')

    fout.close()

if len(sys.argv[:])!=4:
    print "python keepIDRowsByFirstColumn.py json_file(input) output keep_list(first_column_id)"
    sys.exit()

listfile = sys.argv[3]
keep_dic  = listing (listfile)

inputfile = sys.argv[1]
fin = open(inputfile,'r')
obj = json.load(fin)
fin.close()
outputfile = sys.argv[2]

keepProbes (obj, outputfile, keep_dic)
