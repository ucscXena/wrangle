#!/usr/bin/env python

import sys,string,copy
import CGData.RefGene
import CGData.GeneMap
import segToProbeMap

if __name__ == "__main__":
    
    if len(sys.argv[:])!=4:
        print "python mapSegToGeneMatrix.py genomicsSegmentIn regGene GeneLevelMatrixOut\n"
        sys.exit()
    refgene = CGData.RefGene.RefGene()
    refgene.load( sys.argv[2] )

    #* b for cnv 
    probeMapper = CGData.GeneMap.ProbeMapper('b')

    fin =open(sys.argv[1],'r')
    genes= {}
    samples={}
    matrix=[]  #sample then gene
    for gene in refgene.get_gene_list():
        genes[gene]=len(genes)

    Ngene = len(genes.keys())
    oneSample=[]    
    for i in range(0, Ngene):
        oneSample.append([]);

    count =0
    while 1:
        count = count+1
#        print count
        line =string.strip(fin.readline())
        if line =="":
            break
        if line[0]=="#":
            continue
        tmp = string.split(line,"\t")
        if len(tmp)!= 6:
            continue
        seg = segToProbeMap.probeseg("", tmp[1], int(tmp[2]), int(tmp[3]),tmp[4])
        sample = tmp[0]
        value = float(tmp[5])
        if sample not in samples:
            samples[sample]=len(samples)
            matrix.append(copy.deepcopy(oneSample))
            
        hits={}
        for hit in probeMapper.find_overlap( seg, refgene ):
            gene = hit.name
            if gene in hits:
                continue
            hits[gene]=0
            matrix[samples[sample]][genes[gene]].append(value)
    fin.close()

    fout =open(sys.argv[3],'w')
    sample_list =samples.keys()
    fout.write("sample\t"+string.join(sample_list,"\t")+"\n")
    for gene in genes.keys():
        fout.write(gene)
        for sample in sample_list:
            list = matrix[samples[sample]][genes[gene]]
            if len(list)==0:
                average =0
            elif len(list)==1:
                average = list[0]
                average =round(average,3)
            else:
                total=0.0
                for value in list:
                    total = total +value
                average = total/len(list)
                average =round(average,3)
            fout.write("\t"+str(average))
        fout.write("\n")
    fout.close()

