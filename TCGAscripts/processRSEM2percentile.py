import string, os, sys, math

#converting rsem gene or isoform expression results to percentile ranking, rank each gene/transcript by its expression in FPKM/TPM within a single sample, report the percentile ranking

DATA_COL=6 #FPKM

#parse rsem output get FPKM, return dataDic key:FPKM set to NA when there is a duplicated key with different values
def parseInputFPKM(file):
    fin=open(sys.argv[1],'r')
    fin.readline() #header
    dataDic={} # only genes transcripts with one estimates are in the dictionary 
    badIds=[] # genes or transcripts with multiple estimates

    while 1:
        line =fin.readline()
        line =string.strip(line)
        if line =="":
            break
        data = string.split(line,'\t')
        id = data[0]
        value = data[DATA_COL]
        if id not in dataDic:
            try:
                dataDic[id]=float(value)
            except:
                dataDic[id]="NA"
        else:
            try:
                if float(value) == dataDic[id]:
                    pass
                elif id not in badIds:
                    badIds.append(id)
            except:
                dataDic[id]="NA"

    fin.close()
    
    for key in dataDic:
        if key in badIds:
            dataDic[key]="NA"
    
    return dataDic

def percentileRANK (dataDic):
    #sort values within one sample
    values = dataDic.values()
    values.sort()

    #build rank and percentile
    rankDic={}
    c=0
    markV=0
    markR=0

    for v in values:
        if v =="NA":
            break
        if v>0:
            c =c +1
            if markV==0:
                markV=v
            if markR==0:
                markR=c
            if v == markV:
                pass
            else:
                rank = (c +markR)/2.0
                value = markV
                rankDic[value]=rank

                markV=v
                markR=c
        else:
            rankDic[v]=0

    #last one
    rank = (c +markR)/2.0
    value = markV
    rankDic[value]=rank

    #the biggest rank
    N=rank

    dataRankDic={}
    for key in dataDic:
        value = dataDic[key]
        if value in rankDic :
            dataRankDic[key]=rankDic[value]/N*100
        else:
            dataRankDic[key]="NA"
    return dataRankDic

def log2Plus1 (dataDic):
    dataLog2Dic={}
    for key in dataDic:
        value = dataRankDic[key]
        if value != "NA":
            dataLog2Dic[key]= math.log(value+1,2)
        else:
            dataLog2Dic[key]= value
    return dataLog2Dic

def outputDic (outputfile, dataDic,sampleName):
    fout= open(outputfile,'w')
    fout.write("id\t"+sName+"\n")
    for key in dataRankDic:
        fout.write(key+"\t"+str(dataDic[key])+"\n")
    fout.close()

if __name__ == "__main__" :
    if len(sys.argv[:])!= 5:
        print "python processRSEM2percentile.py rsem_result_file percentileFileOut log2(FPKM+1)Fileout sampleName\n"
        sys.exit()

    inputfile = sys.argv[1]
    dataDic = parseInputFPKM(inputfile)
    dataRankDic = percentileRANK (dataDic)
    dataLog2Dic = log2Plus1 (dataDic)

    #sampleName
    sName = sys.argv[4]

    #output percentileRank
    outputfile = sys.argv[2]
    outputDic (outputfile, dataRankDic, sName)

    #output log2(FPKM+1)
    outputfile = sys.argv[3]
    outputDic (outputfile, dataLog2Dic, sName)



