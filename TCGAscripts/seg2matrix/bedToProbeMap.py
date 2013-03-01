#!/usr/bin/env python

import sys,string
import CGData.GenomicSegment
import CGData.SegToMatrix
import CGData.RefGene
import CGData.GeneMap

class beds:    
    def __init__(self):
        self.probes = []
        
    def load (self, handle): #handle bed6
        fin =open(handle,'r')
        while 1:
            line =string.strip(fin.readline())
            string.split(line)
            if line =="":
                break
            tmp = string.split(line,"\t")
            p = probeseg(tmp[3], tmp[0], int(tmp[1]), int(tmp[2]),tmp[5])
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
        print "python bedToProbeMap.py bedInput refGene(eg hg18) probeMapOut\n"
        sys.exit()
        
    probes=beds()
    probes.load(sys.argv[1])
    
    refgene = CGData.RefGene.RefGene()
    refgene.load(sys.argv[2])
    
    #nrefgene = CGData.GeneMap.filter_longest_form(refgene)

    #python pro.py seg refGene_hg18 matrixout
    handle = open(sys.argv[3], "w")
    probeMapper = CGData.GeneMap.ProbeMapper('g')
    for probe in probes.probes:
        hits = []
        for hit in probeMapper.find_overlap( probe, refgene ):
            if hit.name not in hits:
                hits.append(hit.name)
        handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (probe.name, ",".join(hits), probe.chrom, probe.chrom_start, probe.chrom_end, probe.strand))
    handle.close()
    
    
