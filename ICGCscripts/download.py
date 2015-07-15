#!/usr/bin/env python
"""
Download ICGC datasets.
"""

import sys, os, argparse, string, json
from icgcLib import *

def downloadUrl(project, repoName):
    return 'https://dcc.icgc.org/api/v1/download?fn=/' + repInfo['release'] + '/Projects/' + project + '/' + repoName+"."+project+'.tsv.gz'

def downloadIt(projects, dataTypes):
    for p in projects:
        for t in dataTypes:
            url = downloadUrl(p, t)
            file = t + '.' + p + '.tsv.gz'
            try:
                resp = subprocess.check_output(['curl', '--silent', '--fail', '-o', file, url])
            except subprocess.CalledProcessError:
                continue # TODO assuming the file does not exist

            try:
                resp = subprocess.check_output(['curl', '--silent',  '-o',file, url])
            except subprocess.CalledProcessError:
                error('Curl call failed for file: ' + file)
                continue

            info('downloaded: ' + file)

def downloadOriginals():
    """
    Download *.gz files from the repository to the relative directory specified.
    Files are selected based on icgcDataTypes and projects.
    """
    os.chdir(dirOut)

    downloadIt(projects, icgcDataTypes)
    #downloadIt(tcgaProjects, tcgaDataTypes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download ICGC datasets.')
    parser.add_argument('--dirOut', metavar="Directory to write original files from ICGC.", default=dirs.orig)
    args = parser.parse_args()
    dirOut = args.dirOut

    initIcgcLib()
    downloadOriginals()

    os.system("gunzip -f "+dirs.orig+"/*gz")
