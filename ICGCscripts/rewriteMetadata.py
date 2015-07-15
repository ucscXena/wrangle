#!/usr/bin/env python
"""
Rewrite json metadata with different values.
"""
import sys, os, string, csv, json, argparse
from icgcLib import *

def oneFile(fileIn, proj, type):
    if not os.path.exists(fileIn):
        return

    obj={}
    # read in the entire file
    if os.path.exists(fileIn+".json"):
        jsonString = ''
        with open(fileIn +".json", 'rU') as fIn:
            for row in fIn:
                jsonString += row
        obj = json.loads(jsonString)

    obj = buildDatasetCoreMetadata(fileIn, proj)

    with open(fileIn+".json", 'w') as fOut:
        fOut.write(json.dumps(obj, indent=2) + '\n')

def cohortFile(fileIn, proj):
    obj={}
    # read in the entire file
    if os.path.exists(fileIn+".json"):
        jsonString = ''
        with open(fileIn +".json", 'rU') as fIn:
            for row in fIn:
                jsonString += row
        obj = json.loads(jsonString)

    xenaCohort = repInfo['name'] + ' ' + cohortInfo[proj]['xenaSuffix'] +" (" + proj+")"

    obj['description'] = repInfo['cohortDescrPrefix'] + xenaCohort + repInfo['descrSuffix']
    obj['name'] = xenaCohort

    # write out the new json
    with open(fileIn+".json", 'w') as fOut:
        fOut.write(json.dumps(obj, indent=2) + '\n')


def main():
    for p in projects:
        for t in icgcDataTypes:
            fileIn = dirs.xena + t + '.' + p + '.tsv'
            oneFile(fileIn, p, t)
        cohortFile(dirs.xena+'cohort.' + p, p)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rewrite json metadata with different values. No parms.')
    main()
