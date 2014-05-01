#!/usr/bin/env python

import sys,string
import CGData.GenomicSegment
import CGData.SegToMatrix
import CGData.RefGene
import CGData.GeneMap

class segs:    
    def __init__(self):
        self.probes = []
        
    def load (self, handle): #handle bed6
        fin =open(handle,'r')
        while 1:
            line =string.strip(fin.readline())
            if line =="":
                break
            if line[0]=="#":
                continue
            tmp = string.split(line,"\t")
            if len(tmp)== 5:
                p = probeseg(tmp[0], tmp[1], int(tmp[3]), int(tmp[4]),tmp[2])
                self.probes.append(p)
        fin.close()


class probeseg:
    def __init__(self, name, chrom, chrom_start, chrom_end, strand):
        self.name = name
        self.chrom = chrom
        self.chrom_start = chrom_start
        self.chrom_end = chrom_end
        self.strand = strand


if __name__ == "__main__":
    #python pro.py bed refGene_hg18 probeMap
    if len(sys.argv) != 4:
        print "python segToProbeMapGeneExpArray.py segInput(name,chr,strand,start,end) refGene(eg hg18) probeMapOut\n"
        sys.exit()
        
    probes=segs()
    probes.load(sys.argv[1])
    
    refgene = CGData.RefGene.RefGene()
    refgene.load(sys.argv[2])
    
    #nrefgene = CGData.GeneMap.filter_longest_form(refgene)
    #python pro.py seg refGene_hg18 matrixout

    handle = open(sys.argv[3], "w")
    probeMapper = CGData.GeneMap.ProbeMapper('b')
    for probe in probes.probes:
        hits = []
        for hit in probeMapper.find_overlap( probe, refgene ):
            if hit.name not in hits:
                hits.append(hit.name)
        handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (probe.name, ",".join(hits), probe.chrom, probe.chrom_start, probe.chrom_end, probe.strand))
    handle.close()
    
    
