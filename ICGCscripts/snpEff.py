#!/usr/bin/env python
'''
Run snpEff on the simple somatic mutation data. This was the latest version of
snpEff as for 2/28/2015.
'''
import sys, os, string, math, argparse, subprocess
from icgcLib import *

def runSnpEffOnce(inFile):
    infoTime('Processing throught snpEff.')

    # run snpEff
    cmd = 'java -Xmx4g -jar snpEff/snpEff.jar hg19 -c snpEff/snpEff.config -formatEff -quiet -sequenceOntology -noStats -noLog -fileList '
    #cmd = 'java -Xmx4g -jar snpEff/snpEff.jar hg19 -c snpEff/snpEff.config -formatEff -quiet -sequenceOntology -noStats -noLog -no-downstream -no-upstream -no-intergenic -fileList '
    #cmd = 'java -Xmx4g -jar ' + dirs.snpEff + 'snpEff.jar hg19 -c ' + dirs.snpEff + 'snpEff.config -formatEff -quiet -sequenceOntology -noStats -fileList '
    cmd += inFile
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        error('snpEff failed, ' + str(e.returncode) + ': ' + e.output)
        return 1
    infoTime('SnpEff complete.')

def lookup(chr, pos):
    # hg19bases/twoBitToFa hg19bases/hg19.2bit:chr6:31922439-31922452 tmp
    outFile = rootDir + 'hg19bases/tmp'
    cmd = rootDir + 'hg19bases/twoBitToFa ' + rootDir + 'hg19bases/hg19.2bit:chr'
    cmd += chr + ':' + str(pos-1) + '-' + str(pos) + ' ' + outFile
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        error('lookup failed for pos: ' + str(pos) + ', chr: ' + str(chr) + ' with: ' + str(e.returncode) + ': ' + e.output)
        return 1
    with open(outFile,'rU') as fIn:
        r = fIn.next()
        output = fIn.next()
    return string.upper(output[0])

def fixUpInsertionDeletion(r):
    # Handle the case where ref or alt is not given, an insertion or deletion.
    #per https://samtools.github.io/hts-specs/VCFv4.2.pdf, page 4, section 1.4.1.4
    chr = r[0]
    pos = r[1]
    ref = r[3]
    alt = r[4]
    if ref in ['-',".",""]:
        pos = int(pos)
    elif alt in ['-',".",""] and pos not in ["1",1]:
        pos = int(pos) - 1
    else:
        return 1

    val = lookup(chr, pos)

    if val == 1:
        return 1
    if ref in ['-',".",""]:
        ref = val
        alt = val + alt
    elif alt in ['-',".",""]:
        ref = val + ref
        alt = val

    newRow = [chr, pos, r[2], ref, alt]
    return newRow

def writeSampFile(proj, manifest):
    # for each sample ID, run snpEff

    infoTime('Writing input files.')
    snpProjFile = dirs.snp + proj + '.in'
    sort2File = dirs.tmp + 'sort2.' + proj

    # sort by xenaId to group all the same xenaId's together
    os.system('sort -k3 -o ' + sort2File + ' ' + snpProjFile)

    with open(sort2File,'rU') as fIn:
        fIn = csv.reader(fIn, delimiter='\t')
        prevSamp = ''
        for r in fIn:
            if r[0] == 'MT': # some funny chromosome value
                continue

            if r[2] != prevSamp:
                if prevSamp != '':
                    fOut.close()
                    manifest.write(snpSampFile + '\n')
                prevSamp = r[2]
                snpSampFile = dirs.snp + proj + '.'+ prevSamp + '.vcf'
                fOut = open(snpSampFile, 'w')
                fWriter = csv.writer(fOut, delimiter='\t')
                fWriter.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

            if r[3] == '-' or r[4] == '-':
                r = fixUpInsertionDeletion(r)
                if r == 1:
                    return 1

            fWriter.writerow([r[0], r[1], '.', r[3], r[4]])

        try:
            fOut.close()
        except UnboundLocalError:
            pass

        manifest.write(snpSampFile + '\n')

def oneDataset(proj, manifest):
    filebase = 'simple_somatic_mutation.open.' + proj
    fileIn = dirs.orig + filebase + '.tsv'

    cutFile = dirs.tmp + 'cut.' + filebase
    preSortFile = dirs.tmp + 'preSort.' + filebase
    snpProjFile = dirs.snp + proj + '.in'

    if not os.path.exists(fileIn):
        return

    infoTime('SnpEff prep starting for ' + proj + '.')

    xenaIds = loadSampleIdXref(proj)

    # dump unused columns and write to the cut file without headers
    with open(fileIn,'rU') as fIn, open(cutFile,'w') as fOut:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        for r in fIn:
            fOut.writerow([r['chromosome'], r['chromosome_start'], r['icgc_specimen_id'], r['reference_genome_allele'], r['mutated_to_allele']])

    # replace spec ID with xena ID
    with open(cutFile,'rU') as fIn, open(preSortFile,'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        for r in fIn:
            r[2] = xenaIds[r[2]]
            #if r[3] == '-' or r[3] == '': # REF
            #   r[3] = '.'
            #if r[4] == '-' or r[4] == '': # ALT
            #    r[4] = '.'
            fOut.writerow(r)

    # sort to remove duplicates using the key: chr + start + xenaID
    infoTime('Removing duplicates.')
    cmd = 'sort -u -o ' + snpProjFile + ' ' + preSortFile
    os.system(cmd)

    # write one file per sample ID
    retCode = writeSampFile(proj, manifest)
    if retCode == 1:
        error('snpEff failed for ' + proj + '.')
        return 1

    infoTime('snpEff input files written for ' + proj + '.')

def snpEff():
    manifestName = dirs.snp + 'snpEffManifest'
    manifest = open(manifestName, 'w')
    for p in projects:
        retCode = oneDataset(p, manifest)
        if retCode == 1:
            return
    manifest.close()
    runSnpEffOnce(manifestName)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run snpEff on the simple somatic mutation data. No parms.')
    snpEff()
