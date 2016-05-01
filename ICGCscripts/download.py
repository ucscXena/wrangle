#!/usr/bin/env python

bigDir = "/inside/depot/icgcFiles/"
smallDir = "/data/TCGA/icgcFiles/"
release = "release_20"

import os, string, subprocess
import config

def downloadProjectUrl(project, dataType):
    #https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/PACA-CA/copy_number_somatic_mutation.PACA-CA.tsv.gz
    download ="https://dcc.icgc.org/api/v1/download?fn=/"
    return download + release  + '/Projects/' + project + '/' + dataType +"."+project+'.tsv.gz'

def downloadSummaryUrl(dataType):
    #https://dcc.icgc.org/api/v1/download?fn=/release_20/Summary/donor.all_projects.tsv.gz
    download ="https://dcc.icgc.org/api/v1/download?fn=/"
    return download + release  + '/Summary/' + dataType +".all_projects.tsv.gz"

def curlDownload(outdir, file, url):
    try:
        resp = subprocess.check_output(['curl', '--silent', '--fail', '-o', outdir + file, url])
    except subprocess.CalledProcessError:
        return  # TODO assuming the file does not exist
    
    try:
        print url
        print outdir +file
        resp = subprocess.check_output(['curl', '--silent',  '-o',outdir + file, url])
    except subprocess.CalledProcessError:
        error('Curl call failed for file: ' + file)
        return 
    print 'downloaded: ' + file


def downloadOriginals(projects, dataTypes):
    for t in dataTypes:
        for p in projects:
            url = downloadProjectUrl(p, t)
            file = t + '.' + p + '.tsv.gz'
            outdir = smallDir
            if string.find(t,"meth_")!=-1:
                outdir = bigDir            
            curlDownload(outdir, file, url)
        
def downloadSummary(dataTypes):
    for t in dataTypes:
        url = downloadSummaryUrl(t)
        file = t + '.all_projects.tsv.gz'
        outdir = smallDir
        curlDownload(outdir, file, url)

if __name__ == '__main__':
    #downloadOriginals(config.projects, config.icgcDataTypes)
    downloadSummary(config.icgcDataTypes)
    os.system("gunzip -f "+smallDir+"*gz")
    os.system("gunzip -f "+bigDir+"*gz")
