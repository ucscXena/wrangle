#!/usr/bin/env python
"""
Generate the panCan datasets: clinical, RNAseq expression & an abbreviated positional mutation data.
The positional mutation data is only used to create xena IDs.
"""
import sys, os
from icgcLib import *

icgcDataTypes = [ # only the dataset types of interest
    'clinical',
    'exp_seq',
    'simple_somatic_mutation.open',
]

def findFields(fIn, type):
    if type == 'clinical':
        fields = fIn.fieldnames[:] + ['_COHORT']
    elif type == 'exp_seq':
        fields = ['icgc_donor_id', 'normalized_read_count', 'icgc_specimen_id', 'gene_id']
    elif type == 'simple_somatic_mutation.open':
        fields = ['icgc_donor_id', 'icgc_specimen_id'] # only for IDs
    return fields

def oneProjectFile(fileIn, fOut, wroteHead, proj, type):
    infoTime('Processing ' + fileIn)
    with open(fileIn,'rU') as fIn:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fields = findFields(fIn, type)
        if not wroteHead:
            fOut.writerow(fields)
            wroteHead = True
        if type == 'clinical':
            _COHORT = xenaCohortName(proj)
            for r in fIn:
                if r['icgc_donor_id'] != 'icgc_donor_id':
                    r['_COHORT'] = _COHORT
                    fOut.writerow(map(lambda x: r[x], fields))
            return

        for r in fIn:
            if r['icgc_donor_id'] != 'icgc_donor_id':
                fOut.writerow(map(lambda x: r[x], fields))
        return wroteHead

def main():
    infoTime('Running panCanPrep.py')
    for t in icgcDataTypes:
        fileOut = dirs.orig + t + '.ICGC.tsv'
        with open(fileOut, 'w') as fOut:
            fOut = csv.writer(fOut, delimiter='\t')
            wroteHead = False
            for p in projects:
                fileIn = dirs.orig + t + '.' + p + '.tsv'
                if os.path.exists(fileIn) and p != 'ICGC':
                    wroteHead = oneProjectFile(fileIn, fOut, wroteHead, p, t)

if __name__ == '__main__':
    main()


