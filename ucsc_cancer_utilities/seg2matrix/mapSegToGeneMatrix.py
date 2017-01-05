#!/usr/bin/env python

import sys,string,os,copy
import uuid

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/CGData/')

import RefGene
import GeneMap
import segToProbeMap

if __name__ == "__main__":
    if len(sys.argv[:])!=5:
        print "python mapSegToGeneMatrix.py genomicsSegmentIn refGene GeneLevelMatrixOut NORMAL_CNV\n"
        sys.exit()

    refgene = RefGene.RefGene()
    refgene.load( sys.argv[2] )

    NORMAL_CNV=sys.argv[4]

    #* b for cnv
    probeMapper = GeneMap.ProbeMapper('b')

    genes= {}
    samples={}
    for gene in refgene.get_gene_list():
        genes[gene]=len(genes)

    Ngene = len(genes)
    oneSample=[]
    for i in range(0, Ngene):
        oneSample.append('');

    print "genes: ", len(genes)

    # samples
    tmpfile = str(uuid.uuid4())
    os.system("cut -f 1 " + sys.argv[1] + " | sed 1d | sort |uniq > " + tmpfile)
    fin = open(tmpfile, 'r')
    for line in fin.readlines():
        sample = string.strip(line)
        samples[sample]=len(samples)
    fin.close()
    os.system("rm " + tmpfile)

    print "samples: ", len(samples)

    matrix= ([''] * Ngene) * len(samples)  #sample then gene
    matrix_weight =([''] * Ngene) * len(samples)

    fin =open(sys.argv[1],'r')
    count =0
    while 1:
        count = count+1
        line =fin.readline()
        if line =="": # end of file
            break
        if count ==1:
            continue #ignore the first line
        line = string.strip(line)
        if line == "": # empty line
            continue
        if line[0]=="#":
            continue
        tmp = string.split(line,"\t")
        if len(tmp)!= 5:
            continue
        seg = segToProbeMap.probeseg("", tmp[1], int(tmp[2]), int(tmp[3]),".")
        sample = tmp[0]
        value = float(tmp[4])
        hits={}
        for hit in probeMapper.find_overlap( seg, refgene ):
            gene = hit.name
            if gene in hits:
                continue
            hits[gene]=0
            gene_length = hit.chrom_end - hit.chrom_start + 1
            overlap_start = max(int(tmp[2]), hit.chrom_start)
            overlap_end =  min(int(tmp[3]), hit.chrom_end)
            overlap_length = overlap_end - overlap_start + 1
            weight = overlap_length / float(gene_length)

            if matrix[samples[sample] * Ngene + genes[gene]] == '':
                matrix[samples[sample] * Ngene + genes[gene]]  = 0.0
                matrix_weight[samples[sample] * Ngene + genes[gene]] = 0.0
            matrix[samples[sample] * Ngene + genes[gene]] += value * weight
            matrix_weight[samples[sample] * Ngene + genes[gene]] += weight
    fin.close()

    print "segments: ", count

    fout =open(sys.argv[3],'w')
    sample_list =samples.keys()
    fout.write("sample\t"+string.join(sample_list,"\t")+"\n")
    for gene in genes.keys():
        fout.write(gene)
        for sample in sample_list:
            total= matrix[samples[sample] * Ngene + genes[gene]]
            if total == '':
                average = NORMAL_CNV
            else:
                t_weight =matrix_weight[samples[sample] * Ngene + genes[gene]]
                average = total / t_weight
                average =round(average,6)
            fout.write("\t"+str(average))
        fout.write("\n")
    fout.close()

