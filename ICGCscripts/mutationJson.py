'''
Generate json metadata for a mutation file.
'''
import sys, os, string, csv, json
from icgcLib import *

def buildRow(row, meta, stat):
    """
    Gather metadata for this row.
    """
    def value(field):
        if field in row:
            return row[field]
        else:
            return None

    fieldMetaRow(row, meta) # collect metadata before we modify the row

def genomicVector(fileIn, fileOut, fields, stat):
    """
    Process a genomicVector data file.
    """
    with open(fileIn,'rU') as fIn:
        fIn = csv.DictReader(fIn, delimiter='\t')
        row = fIn.next()
        meta = fieldMetaInit(fIn.fieldnames, row, fields)
        outRows = buildRow(row, meta, stat)

        for row in fIn:
            outRows = buildRow(row, meta, stat)

    return meta

def mutationJson(fileIn, fileOut, unused, stat):
    """
    Generate json metadata from an original mutation file.
    """
    fields = ['chromosome', 'chromosome_start', 'icgc_specimen_id', 'reference_genome_allele', 'mutated_to_allele']
    return genomicVector(fileIn, fileOut, fields, stat)
