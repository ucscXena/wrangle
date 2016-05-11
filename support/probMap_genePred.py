import string,sys

#return a dictionary of  hugo:[{(id,)(chr,), (start,), (end,)},]
def parseGenePredToGene(fin):
    chromDic={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        id =data[1]
        chr=data[2]
        strand =data[3]
        start = data[4]
        end =data[5]
        hugo =data[12]
        item ={}
        item['id']=id
        item['chr']=chr
        item['strand']=strand
        item['start']= float(start)
        item['end']=float(end)

        if chr not in chromDic:
            chromDic[chr]={}

        dic = chromDic[chr]

        if hugo!="":
            if hugo not in dic:
                dic[hugo]=[]
            dic[hugo].append(item)

    fin.close()
    return chromDic


def parseProbeMapToGene(fin):
    fin.readline()
    return parseGenePredToGene(fin)

