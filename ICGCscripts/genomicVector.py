"""
Process a genomicVector file. Unused in ICGC r18 wrangling.
"""
import sys, os, string, csv, json
from icgcLib import *

def removeDupPostions(pos):
    delim = '|'
    if string.count(pos, delim) == 0:
        return [pos]
    info('Duplicate position: ' + pos)
    poss = list(set(pos.split(delim))) # make them unique
    poss.sort() # to make test compares easy
    return poss

def buildRow(row, meta, stat):
    """
    Prepare a data row for writing.
    """
    def value(field):
        if field in row:
            return row[field]
        else:
            return None

    fieldMetaRow(row, meta) # collect metadata before we modify the row
    row['chromosome'] = 'chr' + row['chromosome']

    # handle multiple genes: output one row per gene
    hGenes = []
    eGenes = removeDupPostions(row['gene_affected'])
    for g in eGenes:
        hGenes.append(ensemblToHugo(g, stat))
    genes = filter(lambda x: x != None, hGenes)
    i = meta['used'].index('gene_affected')
    rowOut = map(value, meta['used'])

    return map(lambda x: rowOut[:i] + [x] + rowOut[i+1:], genes)

def genomicVector(fileIn, fileOut, headers, fields, stat):
    """
    Process a genomicVector data file.
    """
    sampleId = 'icgc_specimen_id'
    with open(fileIn,'rU') as fIn, open(fileOut, 'w') as fOut:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')

        fOut.writerow(headers) # not yet
        row = fIn.next()
        meta = fieldMetaInit(fIn.fieldnames, row, fields)
        outRows = buildRow(row, meta, stat)
        for r in outRows:
            fOut.writerow(r)

        for row in fIn:
            outRows = buildRow(row, meta, stat)
            for r in outRows:
                fOut.writerow(r)

    return meta

def copy_number(fileIn, fileOut, stat):
    """
    Process a copy number data file.
    """
    headers = ['sample',           'chr',        'start',            'end',            'gene',          'reference', 'alt',               'effect']
    fields  = ['icgc_specimen_id', 'chromosome', 'chromosome_start', 'chromosome_end', 'gene_affected', '',          '',                  'mutation_type']
    return genomicVector(fileIn, fileOut, headers, fields, stat)
