#!/usr/bin/env python
"""
Generate histogram data for a ICGC data file.
"""
import sys, os, argparse, csv, math
from icgcLib import *
from oneHistogram import *

if __name__ == '__main__':
    if len(sys.argv) == 1 or (len(sys.argv) > 1 and sys.argv[1] != 'test'):
        parser = argparse.ArgumentParser(description='Generate histogram data for one ICGC data file.')
        parser.add_argument('fileIn', metavar='Input file name, minus .tsv or tsv.gz extension.')
        parser.add_argument('--dirIn', metavar='Directory of input file.', default=dirs.xena)
        parser.add_argument('--dirOut', metavar='Directory of intermediate file.', default=dirs.hist)
        parser.add_argument('--silent', metavar='Log only to file, not sysout.', default=False)
        parser.add_argument('--all', metavar='Select all samples rather than 10% at random.', default=False)
        args = parser.parse_args()
        makeSilent(args.silent)

        oneHistogram(args.fileIn, args.dirIn, args.dirOut, args.all)
