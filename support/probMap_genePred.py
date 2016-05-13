import string,sys

#return a dictionary of  {chrom:
#                                {hugo:[{(id,),(chr,), (strand,),(start,), (end,)},],}}
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


#id	gene	chrom	chromStart	chromEnd	strand
#ENSG00000223972.5	DDX11L1	chr1	11869	14409	+

def parseProbeMapToGene(fin):
    fin.readline()
    dic={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        id =data[0]
        chr=data[2]
        strand =data[5]
        start = data[3]
        end =data[4]
        hugo =data[1]
        item ={}
        item['id']=id
        item['chr']=chr
        item['strand']=strand
        item['start']= float(start)
        item['end']=float(end)
        item['hugo']= string.split(hugo,','):

        if id not in dic:
            dic[id]=[]
            dic[id]=item
            
    fin.close()
    return dic
