#!/usr/bin/env python
"""
Process a positional xena-ready file to produce a gene-level mutation file.
"""
import sys, os, string, csv, json, subprocess, commands, math
from icgcLib import *
from cfg import *

ones = [ # effects with a value of 1
    'frame',
    'missense',
    'non_coding_exon_variant',
    'splice_acceptor',
    'splice_donor',
    'start',
    'stop_gained',
    'stop_lost',
]

decimalPlaces = 6
igene = 0
isample = 1
ieffect = 2

for i, e in enumerate(ones):
    ones[i] = e.lower().strip()

class Accumulators:
    def __init__(s):
        s.vals = {}
        s.samps = set()
        s.genes = set()
        s.geneNotFound = []
        s.meta = {}

    def printIt(s):
        print 'ACCUM:', s.samps, s.vals

class PreviousRowData:
    def __init__(s):
        s.val = 0
        s.samp = s.gene = ''

    def printIt(s):
        print 'PREV:', s.gene, s.samp, s.val

def updatePositionValue(p, a, newVal):
    """
    Update the value at a sample's position with the previous row's value.
    """
    if p.samp not in a.vals:
        a.vals[p.samp] = 0
    if float(p.val) == 1:
        a.vals[p.samp] = 1
    p.val = newVal

def saveSampleId(p, a, newVal, newSamp):
    """
    Save a sample ID in the samples list.
    """
    updatePositionValue(p, a, newVal)
    if p.samp not in a.samps:
        error('sample not in sample list: ' + p.samp)
    p.samp = newSamp

def findVal(row):
    val = 0
    if 'effect' not in row or row['effect'] == '':
        return val
    effects = string.split(row['effect'], '+')
    for e in effects:
        e = e.strip().lower()
        for one in ones:
            if one in e:
                val = 1
                break
        if val == 1:
            break
    return val

def writeRow(row, accum, fOut):
    accum.samps.add(row['#sample'])
    if 'gene' not in row or row['gene'] == '':
        return
    gene = row['gene']
    val = findVal(row)
    fOut.writerow([gene, row['#sample'], val])

def extractColumns(fileIn, transFile, fields):
    """
    Extract the relevant columns, initialize the accumulators, set the value 
    based on positive effects list & save the column value in the correct order.
    """
    fx = {}
    a = Accumulators()

    with open(fileIn,'rU') as fIn, open(transFile, 'w') as fOut:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')

        fieldsIn = fIn.fieldnames

        # initialize the accumulators and metadata
        r = fIn.next()
        writeRow(r, a, fOut)
        a.meta = fieldMetaInit(fieldsIn, r, fields)
        a.meta['log2'] = False

        # process each row
        for r in fIn:
            writeRow(r, a, fOut)

    return a

def sortByGene(fileIn, sortedFile):
    """
    Sort the input file by the translated gene. We don't need to sort by sample ID
    because sample ID accumulators will handle samples IDs out of order.
    :param fileIn: unsorted data file name
    :param sortedFile: sorted data file name
    :returns: 0 for success, 1 for failure
    """
    infoTime('Starting sort')
    try:
         ret = subprocess.check_output(['sort', '-k', '1,1', '--temporary-directory=data/tmp', '-o', sortedFile, fileIn])
    except subprocess.CalledProcessError:
        error('File sort failed.')
        return 1
    infoTime('Finished sort')

    return 0

def buildOutRow(p, a, newVal, newSamp, newGene):
    """
    Build the row to be written using the accumulated value
    """
    saveSampleId(p, a, newVal, newSamp)

    def finalVal(samp):
        if samp in a.vals:
            return a.vals[samp]
        else:
            return 0
            #return ''

    a.genes.add(p.gene)
    rowOut = [p.gene] + map(finalVal, list(a.samps))
    p.gene = newGene
    a.vals = {}
    return rowOut

def zeroFill(geneFile, dataFileOut, a):
    """
    For each gene not in the original data, create a new row with values of zero for each sample.
    """
    # make lower case: all genes that already have values
    genes = set()
    for gene in a.genes:
        genes.add(gene.lower())

    # step though the ucsc refGene file
    # XXX it would be faster to load the ucsc refGenes into a set() as lower case
    #     outside of this module, as we do for xenaIds
    with open(geneFile,'rU') as fIn, open(dataFileOut, 'a') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        sampCount = len(a.samps)
        for r in fIn:
            if string.lower(r[0]) not in genes:
                fOut.writerow([r[0]] + [0] * sampCount)

def transformToXena(a, sortedFile, fileOut):
    """
    Transform the sorted file to the xena format.
    The input file is sorted by gene with zero or more values per sample-position.
    :param a: accumulators
    :param sortedFile: sorted file name
    :param fileOut: xena-ready file name
    :returns: nothing
    """
    dataFileOut = dirs.tmp + 'data.' + os.path.basename(fileOut)
    geneFile = 'ucscRefGeneHg19.tsv'
    with open(sortedFile,'rU') as fIn, open(dataFileOut, 'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')

        # initialize the previous row values, the accumulators and metadata
        r = fIn.next()
        p = PreviousRowData()
        p.gene = r[igene]
        p.samp = r[isample]
        p.val = r[ieffect]

        for r in fIn:
            val = r[ieffect]
            sample = r[isample]
            if r[igene] == p.gene:
                if sample == p.samp:
                    updatePositionValue(p, a, val)
                else:
                    saveSampleId(p, a, val, sample)
            else:
                rowOut = buildOutRow(p, a, val, sample, r[igene])
                if rowOut != None:
                    fOut.writerow(rowOut)

        # write last data row,
        rowOut = buildOutRow(p, a, 0, '', '')
        if rowOut != None:
            fOut.writerow(rowOut)

    # Create zero filled entries for genes not in original data.
    zeroFill(geneFile, dataFileOut, a)

    # write the header row
    with open(fileOut, 'w') as fOut:
        fOut.write(string.join(['sample'] + list(a.samps), '\t') + '\n')

    # write the data rows and the zero-filled rows
    os.system('cat ' + dataFileOut + ' >> ' + fileOut)
    os.system('rm ' + dataFileOut)

def mutGene(fileIn, fileOut, unused, unused2, unused3):
    infoTime('Starting mutGene.py')
    proj = findCohort(fileIn)
    base = os.path.basename(fileOut)
    fileIn = dirs.xena + 'simple_somatic_mutation.open.' + proj + '.tsv'
    print fileIn, fileOut

    if not os.path.exists(fileIn):
        return 
    initConfig()
    fields = ['gene', '#sample', 'effect']
    transFile = dirs.xLate + base
    sortedFile = dirs.sort + base

    # extract columns of interest
    accum = extractColumns(fileIn, transFile, fields)

    # sort by the translated position
    ret = sortByGene(transFile, sortedFile)
    if ret != 0:
        return None

    # transform this format into xena-ready format
    transformToXena(accum, sortedFile, fileOut)

    return accum.meta

