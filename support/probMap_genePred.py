import string,sys

#return a dictionary of  hugo:[{(id,)(chr,), (start,), (end,)},]
def parseGenePredToGene(fin):
    dic={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        id =data[1]
        chr=data[2]
        strand =data[3]
        start = data[4]
        end =data[5]
        hugo =data[12]
        data['id']=id
        data['chr']=chr
        data['strand']=strand
        data['start']= start
        data['end']=end
        if hugo!="":
            if hugo in dic:
                dic[hugo].append(data)
            else:
                dic[hugo]=[data]
    fin.close()
    return dic 


def parseProbeMapToGene(file):
    dic={}
    fin = open(file,'r')
    fin.readline()
    parseGenePredToGene(fin)
    return dic 

