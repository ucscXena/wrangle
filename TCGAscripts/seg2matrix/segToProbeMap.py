#!/usr/bin/env python

import sys,string
sys.path.insert(0,"../CGDataNew")

import CGData.GenomicSegment
import CGData.SegToMatrix
import CGData.RefGene
import CGData.GeneMap
import optparse

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
    def printUsage():
        print "python segToProbeMap.py segInput(name,chr,strand,start,end) refGene(eg hg18) probeMapOut --mode=cnv|exp\n"

    if len(sys.argv) != 5:
        printUsage()
        sys.exit()

    parser = optparse.OptionParser()
    parser.add_option("--mode", action="store", type="string", dest="mode")
    (options, args) = parser.parse_args()

    if (not options.mode) or (options.mode not in ["cnv","exp"]):
        printUsage()
        sys.exit()

    probes=segs()
    probes.load(sys.argv[1])

    refgene = CGData.RefGene.RefGene()
    refgene.load(sys.argv[2])

    handle = open(sys.argv[3], "w")
    if options.mode=="cnv":
        probeMapper = CGData.GeneMap.ProbeMapper('b')
    if options.mode=="exp":
        probeMapper = CGData.GeneMap.ProbeMapper('g')

    handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("#id", "gene","chrom","chromStart","chromEnd","strand"))
    for probe in probes.probes:
        hits = []
        for hit in probeMapper.find_overlap( probe, refgene ):
            if hit.name not in hits:
                hits.append(hit.name)
        handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (probe.name, ",".join(hits), probe.chrom, probe.chrom_start, probe.chrom_end, probe.strand))
    handle.close()


