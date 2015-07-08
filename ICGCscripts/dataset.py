#!/usr/bin/env python
"""
Process an ICGC data file, excluding mutation data.
"""
import sys, os, argparse, string, json
from icgcLib import *
from oneDataset import *

def dataset(fileIn, dirIn, dirOut):
    initIcgcLib()
    proj = findCohort(fileIn + '.tsv')
    xenaIds = loadSampleIdXref(proj)
    oneDataset(fileIn, dirIn, dirOut, xenaIds)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process one downloaded ICGC data file.')
    parser.add_argument('fileIn', metavar='Original download file base name, in the form: exp_seq.CLLE-ES')
    args = myArgParse(parser)
    dataset(args.fileIn, args.dirIn, args.dirOut)