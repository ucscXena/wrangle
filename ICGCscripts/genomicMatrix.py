"""
Process a genomicMatrix file.
"""
import sys, os, string, csv, json, subprocess, commands, math
from icgcLib import *
from cfg import *

decimalPlaces = 6
xenaSamp = {}

class Accumulators:
    def __init__(s):
        s.vals = {}
        s.samps = []
        s.poss = set()
        s.possNotFound = []
        s.meta = {}

    def printIt(s):
        print json.dumps(s, indent=4, sort_keys=True)

class PreviousRowData:
    def __init__(s):
        s.val = 0
        s.samp = s.pos = ''

    def printIt(s):
        print json.dumps(s, indent=4, sort_keys=True)

def updatePositionCounters(p, a, newVal):
    """
    Update the counters for a position's mean value with the previous row's values.
    """
    if p.samp in a.vals:
        a.vals[p.samp]['count'] += 1
        a.vals[p.samp]['sum'] += p.val
    else:
        a.vals[p.samp] = {'sum': p.val, 'count': 1}
    p.val = newVal

def saveSampleId(p, a, newVal, newSamp):
    """
    Save a sample ID in the samples list.
    """
    updatePositionCounters(p, a, newVal)
    if p.samp not in a.samps:
        a.samps.append(p.samp)
    p.samp = newSamp

def sampleMean(sum, count, a):
    """
    Find the sample value mean.
    """
    if a.meta['log2']:
        val = math.log((float(sum) / count + 1), 2)
    else:
        val = float(sum) / count
    return val

def buildOutRow(p, a, newVal, newSamp, newPos):
    """
    Build the row to be written from the sample means
    """
    saveSampleId(p, a, newVal, newSamp)

    def valMean(samp):
        if samp in a.vals:
            val = sampleMean(a.vals[samp]['sum'], a.vals[samp]['count'], a)
            return val
        else:
            return ''
            
    if p.pos in a.poss and p.pos != '':
        error('Duplicate position data was ignored: ' + str(p.pos))
        # XXX this should not happen because the file should be sorted
        rowOut = None
    else:
        a.poss.add(p.pos)
        rowOut = [p.pos] + map(valMean, a.samps)
    p.pos = newPos
    a.vals = {}
    return rowOut

def extractColumns(fileIn, transFile, fields, stat):
    """
    Extract fields of interest, translate the spec ID to xena ID, and maybe translate to hugo gene.
    :param fileIn: the filename to be read
    :param transFile: the filename to be written
    :param fields: columns to be extracted from fileIn
    :param stat: statistics
    :returns: accum and prev for success, None for failure
    """
    valN = fields[0]
    sampN = fields[1]
    posN = fields[2]
    fx = {}
    a = Accumulators()

    fx['invalidRow'] = "pos == None or row[posN] == '---' or not isNumber(row[valN])"
    fx['getPos'] = 'row[posN]'

    # extract the relevant columns
    with open(fileIn,'rU') as fIn, open(transFile, 'w') as fOut:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')

        fieldsIn = fIn.fieldnames

        # initialize the previous row values, the accumulators and metadata
        for row in fIn:
            pos = eval(fx['getPos'])
            if not eval(fx['invalidRow']):
                break
                
        samp = row[sampN][:]
        if samp[:4] == 'TCGA':
            samp = samp[:15]
        fOut.writerow([pos, samp, row[valN]])

        repoName= findRepoName(fileIn)
        a.meta = fieldMetaInit(fieldsIn, row, fields)
        if 'normalization_algorithm' in row and 'log' not in row['normalization_algorithm'] and repoName in ["exp_array","exp_seq"]: #'RMA' not in row['normalization_algorithm']:
            a.meta['log2'] = True
            info('log2 transform will be applied')
        else:
            a.meta['log2'] = False
            info('log2 transform will NOT be applied')

        # process each row
        for row in fIn:
            pos = eval(fx['getPos'])
            if eval(fx['invalidRow']):
                continue
            samp = row[sampN][:]
            if samp[:4] == 'TCGA':
                samp = samp[:15]
            fOut.writerow([pos, samp, row[valN]])

            fieldMetaRow(row, a.meta)

    return a

def sortByTranslatedPosition(fileIn, sortedFile):
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

def transformToXena(a, sortedFile, fileOut):
    """
    Transform the sorted file to the xena format.
    :param a: accumulators
    :param sortedFile: sorted file name
    :param fileOut: xena-ready file name
    :returns: nothing
    """
    dataFileOut = dirs.tmp + 'data.' + os.path.basename(fileOut)

    with open(sortedFile,'rU') as fIn, open(dataFileOut, 'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')

        # The input file is sorted by genomic position with zero or more values per sample-position.

        # initialize the previous row values, the accumulators and metadata
        posI = 0
        sampI = 1
        valI = 2

        row = fIn.next()
        p = PreviousRowData()
        p.pos = row[posI]
        if row[sampI][:4] == 'TCGA':
            p.samp = row[sampI]
        else:
            p.samp = xenaSamp[row[sampI]]
        p.val = float(row[valI])

        for row in fIn:
            pos = row[posI]
            val = float(row[valI])
            if row[sampI][:4] == 'TCGA':
                samp = row[sampI]
            else:
                samp = xenaSamp[row[sampI]]
            if pos == p.pos:
                if samp == p.samp:
                    updatePositionCounters(p, a, val)
                else:
                    saveSampleId(p, a, val, samp)

            else:
                rowOut = buildOutRow(p, a, val, samp, pos)
                if rowOut != None:
                    fOut.writerow(rowOut)

        # write last data row,
        rowOut = buildOutRow(p, a, 0, '', '')
        if rowOut != None:
            fOut.writerow(rowOut)

    # write the header row
    with open(fileOut, 'w') as fOut:
        fOut.write(string.join(['sample'] + a.samps, '\t') + '\n')

    os.system('cat ' + dataFileOut + ' >> ' + fileOut)
    os.system('rm ' + dataFileOut)

def genomicMatrix(fileIn, fileOut, fields, xenaIds, stat,):
    """
    Process a genomicMatrix data file.
    :param fileIn: input file name as a relative path
    :param fileOut: xena-ready dataset as a relative path
    :param fields: list of field names in input file to build the xena-ready file
    :param xenaIds: icgc_specimen_id to xenaId dictionary
    :param stats: for collecting statistics
    :return: metadata for this dataset
    """
    global xenaSamp
    initConfig()
    xenaSamp = xenaIds
    
    transFile = dirs.xLate + os.path.basename(fileOut)
    sortedFile = dirs.sort + os.path.basename(fileOut)

    # extract columns of interest, and maybe translate to hugo genes
    accum = extractColumns(fileIn, transFile, fields, stat)

    # sort by the translated position
    ret = sortByTranslatedPosition(transFile, sortedFile)
    if ret != 0:
        return None

    # transform this format into xena-ready format
    transformToXena(accum, sortedFile, fileOut)

    return accum.meta

def gene_expression_RNAseq(fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat):
    print fileIn, fileOut
    """
    Process a gene expression RNAseq data file.
    """
    fields = ['normalized_read_count', 'icgc_specimen_id', 'gene_id']
    return genomicMatrix(fileIn, fileOut, fields, xenaIds_specimen, stat)

def gene_expression_array(fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat):
    print fileIn, fileOut
    """
    Process a gene expression Array data file.
    """
    fields = ['normalized_expression_value', 'icgc_specimen_id', 'gene_id']
    return genomicMatrix(fileIn, fileOut, fields, xenaIds_specimen, stat)

def DNA_methylation(fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat):
    print fileIn, fileOut
    """
    Process a gene expression RNAseq data file.
    """
    fields = ['methylation_value', 'icgc_specimen_id', 'probe_id']
    return genomicMatrix(fileIn, fileOut, fields, xenaIds_specimen, stat)

def miRNA_expression_RNAseq(fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat):
    print fileIn, fileOut
    """
    Process an miRNA expression RNAseq data file.
    """
    # for non-TCGA:
    fields = ['normalized_read_count', 'icgc_specimen_id', 'mirna_id']
    return genomicMatrix(fileIn, fileOut, fields, xenaIds_specimen, stat)
