#!/usr/bin/env python
'''
Parse snpEff output for gene and protein.
'''
import sys, os, string, glob, csv, argparse
from icgcLib import *

def oneSamp(sampFile, proj, parsedFile):
    print 'processing', sampFile
    slash = sampFile.index('/')
    dot1 = sampFile.index('.')
    dot2 = sampFile[dot1+1:].index('.') + dot1 + 1
    samp = sampFile[dot1+1:dot2]
    cmd = 'cat ' + sampFile
    cmd += ' | python /inside/home/jzhu/scripts/vcfXenaData/browserDataMelisssa/somaticMutationsForCavm/scripts/parseSnpEffVcf.py '
    cmd += samp + ' ' + parsedFile
    os.system(cmd)

def snpEffParse():
    for p in projects:
        aGlob = dirs.snp + p + '*.eff.vcf'
        files = glob.glob(aGlob)

        if len(files)==0:
            continue

        parsedFile = dirs.xena + 'simple_somatic_mutation.open.' + p + '.tsv'
        with open(parsedFile,'w') as fOut:
            fOut = csv.writer(fOut, delimiter='\t')
            fOut.writerow(['#sample', 'chr', 'start', 'end', 'reference', 'alt', 'gene', 'effect', 'DNA_VAF', 'RNA_VAF', 'amino_acid'])

        for sampFile in files:
            oneSamp(sampFile, p, parsedFile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse snpEff output for gene and protein. No parms.')
    initIcgcLib()
    snpEffParse()
