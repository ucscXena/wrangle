import string, os, sys

def AddByCol(list_vals):
    l1 = map(int,list_vals[0])

    for j in range (1, len(list_vals)):
        l = map(int, list_vals[j])
        for i in range (0, len(l)):
            l1[i] = l1[i] + l[i]
    return l1

def process(inputfile, outputfile):
    dic ={}

    fin = open(inputfile,'r')
    fout = open(outputfile,'w')
    fout.write(fin.readline())

    for line in fin.readlines():
        fout.write(line)
        data = string.split(line[:-1],'\t')
        id = data[0]
        g1, g2 = string.split(id,"->")
        if g1 not in dic:
            dic[g1]= []
        if g2 not in dic:
            dic[g2]= []
        dic[g1].append(data[1:])
        dic[g2].append(data[1:])
    fin.close()

    genes = dic.keys()
    for gene in genes:
        if len(dic[gene]) == 1:
            fout.write(gene+'\t'+string.join(dic[gene][0],'\t')+'\n')
        else:
            total_list = AddByCol(dic[gene])
            fout.write(gene+'\t'+string.join(map(str, total_list),'\t')+'\n')
    fout.close()

if len(sys.argv[:]) != 3:
    print "python parsePKUfusion.py inputfile outputfile"
    sys.exit()

input = sys.argv[1]
output = sys.argv[2]

process(input, output)
