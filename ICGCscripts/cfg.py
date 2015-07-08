"""
Configure for a new Release from ICGC.
"""

import sys, os, gzip, string, csv, json, datetime, subprocess, socket, time
from datetime import date
from itertools import *
from shutil import copyfile

decimalPlaces = 6

rootDir = '/inside/home/jzhu/cgDataJing/ICGCscripts/' # code
dataDir = '/data/TCGA/icgcFiles/' # data


class Dirs: # directory paths
    def __init__(x):

        x.orig = dataDir + 'orig/' # original files downloaded from ICGC
        x.xena = dataDir + 'xena/' # xena-ready files
        x.ids = dataDir + 'ids/' # ID files to map ICGC specimen IDs to xena IDs
        x.origProbeMaps = dataDir + 'origProbeMaps/' # original probeMap input of hsa.gff2.hg19 from mirbase.org
        x.hist = dataDir + 'hist/' # data files to build histograms

        x.tmp = dataDir + 'tmp/' # temporary for intermediate files
        x.origAgain = dataDir + 'tmp/origAgain/' # original files downloaded from ICGC Again to verify by comparing md5's
        x.snp = dataDir + 'tmp/snp/' # snpEff intermediate files
        x.sort = dataDir + 'tmp/sort/' # intermediate sort files
        x.xLate = dataDir + 'tmp/xLate/' # intermediate files for genomic matrix

        x.log = rootDir + 'log/' # log files

dirs = Dirs()
if not os.path.exists(dataDir):
    os.system('mkdir -p ' + dataDir)
iterableDirs = filter(lambda a: not a.startswith('__'), dir(dirs))
for d in iterableDirs:
    di = eval('dirs.' + d)
    os.system('mkdir -p ' + di)

