#!/usr/bin/env python
"""
Generate xena-ready files for matrix data, and the json files for mutations.
"""
import sys, os, argparse, string, glob
from icgcLib import *
from oneDataset import *

def xenaMatrix(dirIn, dirOut):

    infoTime('Started xenaMatrix.py')

    for p in projects:
        xenaIds_specimen = loadSampleIdXref(p)
        xenaIds_donor = loadDonorIdXref(p)
        for t in icgcDataTypes:
            fileBase = t + '.' + p
            if t != 'mutGene':
                origFile = dirIn + '/' + fileBase + '.tsv'
                if not os.path.exists(origFile) and not os.path.exists(origFile + '.gz'):
                    continue
            oneDataset(fileBase, dirIn, dirOut, xenaIds_specimen, xenaIds_donor)

    infoTime('Finished xenaMatrix.py')

if __name__ == '__main__':
    initIcgcLib()
    xenaMatrix(dirs.orig, dirs.xena)
