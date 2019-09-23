import sys,os, string
import urllib2

sys.path.insert(0, os.path.dirname(sys.argv[0]) + "/../support/")
import probMap_genePred

ensemblBasic_url = os.path.dirname(sys.argv[0]) + "/wgEncodeGencodeBasicV19"
gencodeGene_url = os.path.dirname(sys.argv[0]) + "/gencode.v19.annotation.gene.probemap"
upstream = 200

typeRank ={
    "stop_gained":0,
    "frameshift_variant":0,
    "splice_acceptor_variant":0,
    "splice_donor_variant":0,
    "start_lost":1,
    "stop_lost":1,
    "missense_variant":1,
    "inframe_deletion":1,
    "inframe_insertion":1,
    "disruptive_inframe_insertion":1,
    "disruptive_inframe_deletion":1,
    "5_prime_UTR_premature_start_codon_gain_variant":1,
    "splice_region_variant":1,
    "synonymous_variant":2,
    "initiator_codon_variant":2,
    "stop_retained_variant":2,
    "coding_sequence_variant":2,
    "3_prime_UTR_variant":3,
    "5_prime_UTR_variant":3,
    "upstream_gene_variant":4,
    "downstream_gene_variant":4,
    "exon_variant":4,
    "intergenic_regaion":5,
    "intron_variant":5,
    "intragenic_variant":5
}

def enemblIDtoHugo(fin):
    fin.readline()
    dic ={}
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        ensemblG = string.split(data[0],'.')[0]
        hugo = data[1]
        if len(string.split(hugo,','))>1:
            print data
        dic[ensemblG]=hugo
    fin.close()
    return dic

# returen dic key are ensembl Transcripts
def collectT(fin):
    transcripts ={}
    fin.readline()
    for line in fin.readlines():
        ensemblT = line.split('\t')[1]
        transcripts[ensemblT]=0
    fin.close()
    return transcripts

def parseSNV(fin, fout, emsemblGeneDic, upstreamDis):
    fin = open(infile,'r')
    fout = open(outfile,'w')
    dic ={}

    while 1:
        line = fin.readline()
        if line == '':
            break
        data = string.split(line[:-1],'\t')
        sample= data[0]
        chr = data[1]
        start = data[2]
        end = data[3]
        alt = data[5]
        type = data[6]
        enGene = data[8]
        transcript = data[9]

        #only keep upstreamDis variants:
        #if type =="upstream_gene_variant":
        #    if enGene not in emsemblGeneDic:
        #        print enGene, "bad eneGene"
        #        continue
        #    item = emsemblGeneDic[enGene]
            #if item['strand']=="+":
            #    if int(end) < item['start']-upstreamDis:
            #        continue
            #if item['strand']=="-":
            #    if int(start) > item['end']+upstreamDis:
            #        continue

        # record ={}
        # record[transcript]=data
        id = string.join([sample, chr, start, end, alt, enGene],"_")

        if len(dic) != 0  and dic.keys()[0] != id:
            output (fout, dic)
            dic ={}

        if len(dic) == 0:
        #if id not in dic:
            dic[id]={}

        if transcript not in dic[id]:
            dic[id][transcript]=data
        elif dic[id][transcript][6] != data[6]:
            print id, "old", dic[id][transcript][6], "new", data

    fin.close()
    fout.close()

def output(fout, data):
    for key in data:
        transcripts = data[key].keys()

        #use ensembl canonical
        #foundC =0
        #for t in transcripts:
        #    if t in cTranscripts:
        #        record = data[key][t]
        #        hugo= findHugo(idToHugo, record[8], record[9])
        #        record[8] = hugo
        #        fout.write(string.join(record[:9], "\t")+"\n")
        #        foundC= foundC +1
        #if foundC:
        #    continue
        
        if len(data[key])==1:
            t = data[key].keys()[0]
            record = data[key][t]
            hugo= findHugo(idToHugo, record[8], record[9])
            record[8] = hugo
            fout.write(string.join(record[:9], "\t")+"\n")
            continue

        uniq_v =[]
        for value in data[key].values():
            if value[:9] not in uniq_v:
                uniq_v.append(value[:9])

        if len(uniq_v)==1:
            record = uniq_v[0]
            hugo= findHugo(idToHugo, record[8], data[key].values()[0][9])
            record[8] = hugo
            fout.write(string.join(record, "\t")+"\n")
            continue

        #use ensembl basic
        foundC =0
        record =[]
        for t in transcripts:
            if t in basicTranscripts:
                record = data[key][t]
                foundC= foundC +1

        if foundC==1:
            hugo= findHugo(idToHugo, record[8], record[9])
            record[8] = hugo
            fout.write(string.join(record[:9], "\t")+"\n")
            continue

        #selet best type of mutation within ensembl basic
        minRank = 6
        bestRecord=[]
        for t in transcripts:
            if t in basicTranscripts:
                record = data[key][t]
                type = record[6]
                if type in typeRank:
                    rank = typeRank[type]
                else:
                    print "missing type:", type
                    continue
                if rank < minRank:
                    minRank = rank
                    bestRecord = record

        if bestRecord !=[]:
            hugo= findHugo(idToHugo, bestRecord[8], bestRecord[9])
            record[8] = hugo
            bestRecord[8] = hugo
            fout.write(string.join(bestRecord[:9], "\t")+"\n")
            continue

        #only hits in non-basic transcript
        fout.write(string.join(data[key].values()[0][:6], "\t")+"\t\t\t\n")

def findHugo(idToHugo, gene, transcript):
    if gene in idToHugo:
        return idToHugo[gene]
    if transcript in idToHugo:
        return idToHugo[transcript]
    return ""

if len(sys.argv[:])!= 3:
    print "python icgcSNVcleanup.py input output"
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

fin = open(ensemblBasic_url, 'r')
basicTranscripts = collectT(fin)

fin = open(gencodeGene_url, 'r')
idToHugo = enemblIDtoHugo(fin)

fin = open(gencodeGene_url, 'r')
g = lambda key: string.split(key,'.')[0]
emsemblGeneDic = probMap_genePred.parseProbeMapToGene(fin, g)

####### parse
data = parseSNV(infile, outfile, emsemblGeneDic, upstream)


