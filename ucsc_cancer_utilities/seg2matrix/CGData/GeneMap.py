#!/usr/bin/env python

import csv
import sys
import re

import CGData.ProbeMap
import CGData.GenomicMatrix

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



def genomicSegment2Matrix(genomicSegment, refGene, probeMapper):
    """
    Take a genomicSegment map, compare it against a refGene table,
    and contruct a genomicMatrix
    """
    out = CGData.GenomicMatrix.GenomicMatrix()
    out.init_blank( rows=refGene.get_gene_list(), cols=genomicSegment.get_key_list() )
    for id in genomicSegment.get_key_list():
        for segment in genomicSegment.get_by(id):
            for hit in probeMapper.find_overlap( segment, refGene ):
                out.set_val(row_name=hit.name, col_name=segment.id, value=segment.value )
    return out


def filter_longest_form(refgene):
    """
    take a refgene table and filter multiple gene isoforms down to the longest
    """
    ng = CGData.RefGene.RefGene()
    for g in refgene.get_gene_list():
        longest = None
        length = 0
        for elem in refgene.get_gene(g):
            newLength = elem.chrom_end - elem.chrom_start
            if newLength > length:
                length = newLength
                longest = elem
        ng.add(longest)
    ng.loaded = True
    return ng


def genomicSegment2MatrixNorm(genomicSegment, refGene, probeMapper):
    """
    Given 
    """
    ng = filter_longest_form(refGene)
    #enumerate the col order of the sample ids
    idList = genomicSegment.get_key_list()
    
    geneList = ng.get_gene_list()
    
    tmp = CGData.GenomicMatrix.GenomicMatrix()
    tmp.init_blank( rows=geneList, cols=idList )
    geneHits = {}
    #read through the segment one sample id at a time
    for id in idList:   
        segmentMap = {}
        for segment in genomicSegment.get_by(id):
            for hit in probeMapper.find_overlap( segment, ng ):
                span = float(min(segment.chrom_end, hit.chrom_end) - max(segment.chrom_start, hit.chrom_start)) / float(hit.chrom_end - hit.chrom_start)
                #if hit.name not in segmentMap:
                #    segmentMap[hit.name] = []
                try:
                    segmentMap[hit.name].append(
                        ( span, segment.value )
                    )
                except KeyError:
                    segmentMap[hit.name] = [
                        ( span, segment.value )
                    ]
        
        for gene in segmentMap:
            geneHits[gene] = True
            mapInfo = segmentMap[gene]
            coverage = sum( i[0] for i in mapInfo )
            assert coverage <= 1.0
            value = sum( i[0]*i[1] for i in mapInfo )
            #print coverage, value, value/coverage, segmentMap[gene]
            tmp.set_val(row_name=gene, col_name=id, value=value/coverage)
    
    #now remove the blanks
    out = CGData.GenomicMatrix.GenomicMatrix()
    out.init_blank( rows=geneHits, cols=idList )
    for gene in geneHits:
        for sample in idList:
            out.set_val( row_name=gene, col_name=sample, value=tmp.get_val(row_name=gene, col_name=sample))
    return out



def aliasRemap(genomicMatrix, aliasMap):
    """
    Given a genomicMatrix and an alias map, create a new genomic matrix 
    with the probes from the original matrix remapped to the connected aliases
    from the map
    """
    
    am = {}
    for probe in aliasMap.get_key_list():
        for alias in aliasMap.get_by(probe):
            if alias not in am:
                am[alias.alias] = {}
            am[alias.alias][probe] = True
        
    out = CGData.GenomicMatrix.GenomicMatrix()
    out.init_blank( rows=am.keys(), cols=genomicMatrix.get_col_list() )
    probeMap = genomicMatrix.get_row_map()
    for a in am:
        for sample in genomicMatrix.get_col_list():
            o = []
            for p in am[a]:
                if p in probeMap:
                    o.append( genomicMatrix.get_val( col_name=sample, row_name=p) )
            if len(o):
                out.set_val(col_name=sample, row_name=a, value=sum(o) / float(len(o)))
    
    return out
    
def refGeneLink2ProbeLoc(aliasMap, refGene):
    """
    given an alias map, and a refGene produce a probeMap by connecting 
    alias symbols. Returns the coordinates of the longest form
    """
    out = CGData.ProbeLoc.ProbeLoc()
    out.init_blank()
    
    for probe in aliasMap.get_key_list():
        for link in aliasMap.get_by(probe):
            probe = link.probe
            geneName = link.alias
            sGene = None
            try:
                for gene in refGene.get_gene(geneName):
                    if sGene is None or gene.chrom_end - gene.chrom_start > sGene.chrom_end - sGene.chrom_start:
                        sGene = gene
            except KeyError:
                pass
        if sGene is not None:
            out.insert(probe, { 'probe' : probe, 'chrom' : sGene.chrom, 'strand' : sGene.strand, 'chrom_start' : sGene.chrom_start, 'chrom_end' : sGene.chrom_end })
    return out


    
    
    
