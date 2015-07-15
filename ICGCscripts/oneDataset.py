"""
Inner routines to process a downloaded file, excluding mutation data.
"""
import sys, os, argparse, string, json
from icgcLib import *
from clinicalMatrix import *
from genomicVector import copy_number
from mutationJson import mutationJson
from genomicMatrix import gene_expression_array, gene_expression_RNAseq, DNA_methylation, miRNA_expression_RNAseq
from mutGene import mutGene
from download import *

def oneDataset(fileIn, dirIn, dirOut,  xenaIds_specimen, xenaIds_donor):
    if fileIn[-4:] == '.tsv' or fileIn[-7:] == '.tsv.gz' or fileIn.count('/') > 0:
        error("Input file must be without a directory and without a '.tsv' or '.tsv.gz' extension.")
        return

    fileOut = dirOut + fileIn + '.tsv'
    infoTime('*** processing file: ' + fileIn)

    fileIn = dirIn + fileIn + '.tsv'

    if fileIn.count('mutGene') == 0:
        if not os.path.exists(fileIn):
            if not os.path.exists(fileIn + '.gz'):
                error('Neither of these input files exist: ' + fileIn + ', ' + fileIn + '.gz')
                return

            infoTime('Starting gunzip.')
            try:
                subprocess.check_output(['gunzip', fileIn + '.gz'])
            except subprocess.CalledProcessError:
                error('Gunzip failed.')
                return

            infoTime('Finished gunzip.')
        if not os.path.exists(fileIn):
            return 

    stat = {
        'geneCache': {},
        'emptyPositionCount': 0,
        'emptyGeneFound': False,
        'noHugoGenes': [],
        'totalEntries': 0,
        'geneXlateFx': 'snpEff(eGene)',
    }


    # write the xena-ready dataset file
    dataSubType = findDataSubType(fileIn)

    if dataSubType == 'invalid':
        infoTime('Failed processing file: ' + os.path.basename(fileIn))
        return

    fields = eval(dataSubType['dsFx'] + '(fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat)')
    if fields == None:
        infoTime('Failed processing file: ' + os.path.basename(fileIn))
        return

    cohort = findCohort(fileIn)
    info('writing cohort json for ' + cohort)
    writeCohortMetadata(cohort, dirOut)
    repName= findRepoName (fileIn)
    url = downloadUrlFromFile(cohort, os.path.basename(fileIn))
    LOG = False
    if fields.has_key('meta') and fields.meta.has_key('log2') and fields.meta.log2:
        LOG=True
    fields  = buildDatasetCoreMetadata(repName, url, cohort, LOG)
    writeMetadata(fields, fileOut)
    if 'geneTrans' in fields and fields['geneTrans']:
        uniqueGenes = len(stat['geneCache'])
        info(str(stat['emptyPositionCount']) + ' empty position fields.')
        info(str(stat['totalEntries']) + ' lines written to xena-ready dataset.')
        if uniqueGenes:
            noHugoUnique = sorted(list(set(stat['noHugoGenes'])))
            noHugoUniqueLen = len(noHugoUnique)
            info(str(len(stat['noHugoGenes'])) + ' lines written to xena-ready dataset with no Hugo translation.')
            info(str(uniqueGenes) + ' unique ensembl genes.')
            info(str(noHugoUniqueLen) + ' unique ensembl genes with no Hugo translation')
            info(str(int(round(((float(uniqueGenes) - noHugoUniqueLen) / uniqueGenes) * 100))) + '% gene translation coverage.')
        
    infoTime('Completed processing file: ' + os.path.basename(fileIn) + ', if successful, generated ' + fileOut)
