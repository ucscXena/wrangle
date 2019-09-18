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
    try:
        print url
        print outdir +file
        subprocess.check_output(['wget', '-O', outdir + file, url, '--no-check-certificate'])
        print 'downloaded: ' + file
    except subprocess.CalledProcessError:
        print 'Curl call failed for file: ' + file
        os.system("rm -f "+ outdir + file)

def downloadOriginals(projects, dataTypes):
    print projects
    print dataTypes
    for t in dataTypes:
        for p in projects:
            #if p == "BRCA-EU":
            #    continue
            url = downloadProjectUrl(p, t)
            file = t + '.' + p + '.tsv.gz'
            outdir = config.downloadDir
            try:
                curlDownload(outdir, file, url)
            except:
                print t, p
                pass
def downloadSummary(dataTypes):
    for t in dataTypes:
        url = downloadSummaryUrl(t)
        file = t + '.all_projects.tsv.gz'
        print file
        outdir = config.downloadSummaryDir
        curlDownload(outdir, file, url)

if __name__ == '__main__':
    print config.getProjects()
    downloadOriginals(config.getProjects(), config.icgcDataTypes)
    ##downloadOriginals(["HNSC-US","LUAD-US","BRCA-US"], config.icgcDataTypes)
    downloadSummary(config.icgcDataTypes)
    os.system("gunzip -f "+config.downloadDir+"*gz")
