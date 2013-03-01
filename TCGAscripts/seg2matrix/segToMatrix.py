#!/usr/bin/env python

import sys
import CGData.GenomicSegment
import CGData.SegToMatrix
import CGData.RefGene
import CGData.GeneMap

class matrix_write:    
    def __init__(self, handle):
        self.buff = ""
        self.handle = handle
        self.probes = []
    def write(self, s):
        self.buff += s        
        if s.endswith("\n"):
            tmp = self.buff.split("\t")
            if tmp[0] != "probe":
                tmp2 = tmp[0].split("_")
                p = probeseg(tmp[0], tmp2[0], int(tmp2[1]), int(tmp2[2]))
                self.probes.append(p)
                
            self.handle.write(self.buff)
            self.buff = ""

class probeseg:
    def __init__(self, name, chrom, chrom_start, chrom_end):
        self.name = name
        self.chrom = chrom
        self.chrom_start = chrom_start
        self.chrom_end = chrom_end
        self.strand = "."


if __name__ == "__main__":
    seg = CGData.GenomicSegment.GenomicSegment()
    seg.load(sys.argv[1])
    
    refgene = CGData.RefGene.RefGene()
    refgene.load(sys.argv[2])
    
    #nrefgene = CGData.GeneMap.filter_longest_form(refgene)

    outbase = "out"
    
    if len(sys.argv) > 3:
        outbase = sys.argv[3]
        
    handle = open(outbase + ".matrix", "w")
    m = matrix_write(handle)    
    CGData.SegToMatrix.seg_to_matrix(seg, m)
    handle.close()
    
    handle = open(outbase + ".probeMap", "w")
    probeMapper = CGData.GeneMap.ProbeMapper('b')
    for probe in m.probes:
        hits = []
        for hit in probeMapper.find_overlap( probe, refgene ):
            if hit.name not in hit:
                hits.append(hit.name)
        handle.write("%s\t%s\t%s\t%s\t%s\t.\n" % (probe.name, ",".join(hits), probe.chrom, probe.chrom_start, probe.chrom_end))
    handle.close()
    
    
