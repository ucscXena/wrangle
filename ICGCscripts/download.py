#!/usr/bin/env python

bigDir = "/inside/depot/icgcFiles/"
smallDir = "/data/TCGA/icgcFiles/"

import os, string, subprocess
from icgcLib import projects, icgcDataTypes

def downloadUrl(project, dataType):
    #https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/PACA-CA/copy_number_somatic_mutation.PACA-CA.tsv.gz
    download ="https://dcc.icgc.org/api/v1/download?fn=/"
    release = "release_20"
    return download + release  + '/Projects/' + project + '/' + dataType +"."+project+'.tsv.gz'

def downloadOriginals(projects, dataTypes):
    for p in projects:
        for t in dataTypes:
            url = downloadUrl(p, t)
            file = t + '.' + p + '.tsv.gz'
            outdir = smallDir
            if string.find(t,"meth_")!=-1:
                outdir = bigDir
            try:
                resp = subprocess.check_output(['curl', '--silent', '--fail', '-o', outdir + file, url])
            except subprocess.CalledProcessError:
                continue # TODO assuming the file does not exist

            try:
                print url
                print outdir +file
                resp = subprocess.check_output(['curl', '--silent',  '-o',outdir + file, url])
            except subprocess.CalledProcessError:
                error('Curl call failed for file: ' + file)
                continue

            print 'downloaded: ' + file

if __name__ == '__main__':
    downloadOriginals(projects, icgcDataTypes)
    os.system("gunzip -f "+smallDir+"*gz")
    os.system("gunzip -f "+bigDir+"*gz")
