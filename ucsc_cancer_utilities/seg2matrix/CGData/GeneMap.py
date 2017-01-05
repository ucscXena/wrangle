#!/usr/bin/env python

import csv
import sys
import re

class ProbeMapper(object):
    """
    Class to map the probes. Expects handle to the refGene_hg18.table file
    """

    def __init__(self, mode='g'):
        self.cmp_func = optionMap[mode]

    def find_overlap(self, segment, ref_gene, cmp_func=None):
        """
        Function to find overlaps for a given probe description.
        the cmp_func arg is a function that returns a 'True' or 'False' for
        a given probe description and a gene, examples include 'gene_overlap'
        and 'gene_simple_meth_overlap'
        """
        if cmp_func is None:
            cmp_func = self.cmp_func

        if not ref_gene.has_chrom(segment.chrom):
            return []
        chromList = ref_gene.get_chrom(segment.chrom)

        out = []
        for gene in chromList:
            if cmp_func(segment.chrom_start,
                        segment.chrom_end, segment.strand, gene):
                out.append(gene)
        return out


#
# The set of functions that can be used to do comparisons
#


def block_both_strand(start, end, strand, gene):
    """
    Check is segment is between gene start and end, either strand
    
    **Code 'b'**
    """
    if gene.chrom_end >= start and gene.chrom_start <= end:
        return True
    return False

def block_same_strand(start, end, strand, gene):
    """
    Check is segment is on same strand, between gene start and end
    
    **Code 's'**
    """
    if gene.chrom_end >= start and gene.chrom_start <= end and strand == gene.strand:
        return True
    return False


def exon_same_strand(start, end, strand, gene):
    """
    Check is segment is on same strand, and occurs on an exon
    
    **Code 'g'**
    """

    if gene.strand != strand:
        return False
    for i in range(int(gene.ex_count)):
        if gene.ex_end[i] >= start and gene.ex_start[i] <= end:
            return True
    return False


def exon_both_strand(start, end, strand, gene):
    """
    Check is segment occurs on an exon, on either stand
    
    **Code 'e'**
    """
    for i in range(int(gene.ex_count)):
        if gene.ex_end[i] >= start and gene.ex_start[i] <= end:
            return True
    return False

def block_same_strand_coverage75(start, end, strand, gene):
    if gene.chrom_end >= start and gene.chrom_start <= end and strand == gene.strand:
        cov = min(gene.chrom_end,end) - max(gene.chrom_start, start)
        if float(cov) / float(end - start) > 0.75:
            return True
    return False


def exon_same_strand_coverage75(start, end, strand, gene):
    if strand != gene.strand:
        return False
    cov = 0
    for i in range(int(gene.ex_count)):
        if gene.ex_end[i] >= start and gene.ex_start[i] <= end:
            cov += min(gene.ex_end[i],end) - max(gene.ex_start[i], start)
    if float(cov) / float(end - start) > 0.75:
        return True
    return False


class Intersector:
    def __init__(self, same_strand=True, exons=False, coverage=None, start_rel_cdsStart=None, end_rel_cdsStart=None, start_rel_tss=None, end_rel_tss=None):
        self.same_strand = same_strand
        self.exons = exons
        self.coverage = coverage
        self.start_rel_tss = start_rel_tss
        self.end_rel_tss = end_rel_tss
        self.start_rel_cdsStart = start_rel_cdsStart
        self.end_rel_cdsStart = end_rel_cdsStart
    
    def hit(self, start, end, strand, gene):
        if self.same_strand and strand != gene.strand:
            return False

        if self.exons:
            cov = 0
            for i in range(int(gene.ex_count)):
                if gene.ex_end[i] >= start and gene.ex_start[i] <= end:
                    cov += min(gene.ex_end[i],end) - max(gene.ex_start[i], start)
            if self.coverage is None and cov > 0:
                return True
            else:
                if float(cov) / float(end - start) > self.coverage:
                    return True
        else:
                wStart = gene.chrom_start
                wEnd   = gene.chrom_end
                
                if self.start_rel_cdsStart is not None:                
                    if gene.strand == "+":
                        wStart = gene.cds_start + self.start_rel_cdsStart
                    if gene.strand == "-":
                        wStart = gene.cds_end - self.start_rel_cdsStart
                if self.end_rel_cdsStart is not None:                
                    if gene.strand == "+":
                        wEnd = gene.cds_start + self.end_rel_cdsStart
                    if gene.strand == "-":
                        wEnd = gene.cds_end - self.end_rel_cdsStart
                
                if self.start_rel_tss is not None:
                    if gene.strand == "+":
                        wStart = gene.chrom_start + self.start_rel_tss
                    if gene.strand == "-":
                        wStart = gene.chrom_end - self.start_rel_tss
                if self.end_rel_tss is not None:
                    if gene.strand == "+":
                        wEnd = gene.chrom_start + self.end_rel_tss
                    if gene.strand == "-":
                        wEnd = gene.chrom_end - self.end_rel_tss
                
                cstart = min(wEnd, wStart)
                cend = max(wEnd, wStart)
                if cend >= start and cstart <= end:
                    cov = min(gene.chrom_end,end) - max(gene.chrom_start, start)
                    if self.coverage is None and cov > 0:
                        return True
                    else:
                        if float(cov) / float(end - start) > 0.75:
                            return True
    
        return False

        

###ADD MORE FUNCTIONS HERE


####

###To add options to the command line, map the option character to a function
###for example '-m' maps to gene_simple_meth_overlap

optionMap = {
    "b": block_both_strand,
    "s": block_same_strand,
    "g": exon_same_strand,
    "e": exon_both_strand
}


