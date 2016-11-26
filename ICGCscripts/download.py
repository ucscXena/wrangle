#!/usr/bin/env python

import os, string, subprocess
import config

def downloadProjectUrl(project, dataType):
    download ="https://dcc.icgc.org/api/v1/download?fn=/"
    return download + config.release  + '/Projects/' + project + '/' + dataType +"."+project+'.tsv.gz'

def downloadSummaryUrl(dataType):
    download ="https://dcc.icgc.org/api/v1/download?fn=/"
    return download + config.release  + '/Summary/' + dataType +".all_projects.tsv.gz"

def curlDownload(outdir, file, url):
    #try:
    #    resp = subprocess.check_output(['wget', '-O', outdir + file, url,'--no-check-certificate'])
    #except subprocess.CalledProcessError:
    #    return  # TODO assuming the file does not exist
    
    try:
        print url
        print outdir +file
        resp = subprocess.check_output(['wget', '-O', outdir + file, url, '--no-check-certificate'])
        print 'downloaded: ' + file
    except subprocess.CalledProcessError:
        print 'Curl call failed for file: ' + file
        os.system("rm -f "+ outdir + file)


def downloadOriginals(projects, dataTypes):
    for t in dataTypes:
        for p in projects:
            url = downloadProjectUrl(p, t)
            file = t + '.' + p + '.tsv.gz'
            outdir = config.smallDir
            if string.find(t,"meth_")!=-1:
                outdir = config.bigDir
            curlDownload(outdir, file, url)
        
def downloadSummary(dataTypes):
    for t in dataTypes:
        url = downloadSummaryUrl(t)
        file = t + '.all_projects.tsv.gz'
        print file
        outdir = config.smallDir
        curlDownload(outdir, file, url)

if __name__ == '__main__':
    downloadOriginals(config.getProjects(), config.icgcDataTypes)
    #downloadOriginals(["PACA-AU"], config.icgcDataTypes)
    #downloadSummary(config.icgcDataTypes)
    os.system("gunzip -f "+config.smallDir+"*gz")
    os.system("gunzip -f "+config.bigDir+"*gz")
