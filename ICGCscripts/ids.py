#!/usr/bin/env python
"""
Generate a dictionary with entries of the form:  icgc_specimen_id: xenaId
 input file: orig/*
output file: ids/xref.<proj>
"""
import sys, os, string, math, argparse
from icgcLib import *

cuts = { # columns of interest to cut from original files
    'specimen': '1,5,7,8',
    #'clinical': '1,18,20,25',
    'copy_number_somatic_mutation': '1,3',
    'exp_array': '1,3',
    'exp_seq': '1,3',
    'meth_array': '1,3',
    'mirna_seq': '1,3',
    'simple_somatic_mutation.open': '2,4',
}

def oneDataset(type, proj, catFile):
    fileBase = type + '.' + proj
    fileIn = dirs.orig + fileBase + '.tsv'
    cutFile = dirs.tmp + 'cut.' + fileBase
    sortTypeFile = dirs.tmp + 'sort.' + fileBase
    typeFile = dirs.tmp + 'type.' + fileBase

    if not os.path.exists(fileIn):
        return

    # Extract fields containing only the IDs plus for clinical, also extract
    # icgc_specimen_type and icgc_specimen_processing. Then sort the data.
    os.system('cut -f ' + cuts[type] + ' ' + fileIn + ' > ' + cutFile)
    os.system('sort -u -o ' + sortTypeFile + ' ' + cutFile) # sort on donor, then specId

    # Rewrite the data, inserting the dataType on each entry, and cat the
    # results to a project accumulation file. NOTE: this could be combined with
    # the above.
    with open(sortTypeFile,'rU') as fIn, open(typeFile, 'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        next(fIn)# skip header row
        for row in fIn:
            if type == 'specimen':
                rowOut = [row[1], row[0], type, row[2], row[3]]
            else:
                rowOut = [row[0], row[1], type]

            fOut.writerow(rowOut)
    cmd = 'cat ' + typeFile + ' >> ' + catFile
    os.system(cmd)

def buildSampleId(donor, specType):
    keep = string.letters + string.digits + '_'
    type = ''
    for letter in specType:
        if letter in keep:
            type += letter
        else:
            type += '_'
    type = type.replace('___', '_')
    type = type.replace('__', '_')
    if type[len(type) - 1] == '_':
        type = type[:len(type) - 1]
    return donor + '_' + type

def buildDictionary(proj):
    """
    Generate a dictionary keyed on icgc_specimen_id, referencing its xena ID
    """
    fileIn = dirs.tmp + 'uniq1.' + proj
    orderFile = dirs.tmp + 'order.' + proj
    sortFile = dirs.tmp + 'sort.' + proj
    xrefFile = dirs.ids + 'xref.' + proj
    if not os.path.exists(fileIn):
        print 'no file at:', fileIn
        return

    # put row values in the right order, then sort the file
    with open(fileIn,'rU') as fIn, open(orderFile,'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        for row in fIn:
            # write according to new sort order: donor, specType, dataType, specID
            fOut.writerow([row[0], row[3], row[2], row[1]])

    os.system('sort -u -o ' + sortFile + ' ' + orderFile)

    # Generate a cross-reference file keyed on icgc_specimen_id, referencing our
    # newly constructed sample ID.
    donor = 0
    specType = 1
    dataType = 2
    specId = 3

    with open(sortFile,'rU') as fIn, open(xrefFile,'w') as fOut:
        fIn = csv.reader(fIn, delimiter='\t')
        fOut = csv.writer(fOut, delimiter='\t')
        try:
            r = fIn.next()
            prevSpec = r[specId]
            prevSample = buildSampleId(r[donor], r[specType])
        except StopIteration:
            return
        for r in fIn:
            if r[specId] == prevSpec:
                continue

            # write out the reference row for the previous specimen ID
            fOut.writerow([prevSpec, prevSample])

            prevSpec = r[specId]
            prevSample = buildSampleId(r[donor], r[specType])

        # write out the last reference row
        fOut.writerow([prevSpec, prevSample])

def ids():
    for p in projects:
        catFile = dirs.tmp + 'cat.' + p
        sortFile = dirs.tmp + 'sort.' + p
        os.system('rm -f ' + catFile)
        for t in icgcDataTypes:
            if string.find(t,"donor")==0:
                continue
            oneDataset(t, p, catFile)

        # Sort the cat file containing IDs from all dataTypes for a project
        cmd = 'sort -o ' + sortFile + ' ' + catFile # sort by donor, then specId, then dataType
        os.system(cmd)

        # Add the specimen type and processing found in the clinical, to the
        # other dataType ID entries.
        specTypeFile = dirs.tmp + 'uniq.' + p

        with open(sortFile,'rU') as fIn, open(specTypeFile, 'w') as fOut:
            fIn = csv.reader(fIn, delimiter='\t')
            fOut = csv.writer(fOut, delimiter='\t')
            for row in fIn:
                rowOut = [row[0], row[1], row[2]]
                if row[2] == 'specimen': # clinical is always first in the sort
                    specType = row[3]
                    specProc = row[4]
                fOut.writerow(rowOut + [specType, specProc])

        # Remove clinical entries with no genomic data
        uniqFile = dirs.tmp + 'uniq.' + p
        uniq1File = dirs.tmp + 'uniq1.' + p
        with open(uniqFile,'rU') as fIn, open(uniq1File, 'w') as fOut:
            fIn = csv.reader(fIn, delimiter='\t')
            fOut = csv.writer(fOut, delimiter='\t')
            r = fIn.next()
            prevSpecId = r[1]
            prevDataType = r[2]
            prevRow = r
            for r in fIn:
                if r[1] == prevSpecId:
                    fOut.writerow(prevRow)
                elif prevDataType != 'specimen':
                    fOut.writerow(prevRow)
                prevSpecId = r[1]
                prevDataType = r[2]
                prevRow = r
            if prevDataType != 'specimen':
                fOut.writerow(prevRow)

        buildDictionary(p)

import re
newCUTs= '1,5,7'

def newIDs():
    for proj in projects:
        fileIn = dirs.orig + "specimen." + proj+'.tsv'
        tmpfile = ".tmp"
        fileOut = dirs.ids + "mapping." + proj+'.tsv'
        os.system("cut -f "+ newCUTs + " " +fileIn+" > "+ tmpfile)
        
        fin =open(tmpfile ,'U')
        fout = open(fileOut,'w')
        fout.write(fin.readline())
        
        for line in fin.readlines():
            sID, dID, sType = string.split(line[:-1],'\t')
            newType = re.sub('[^0-9a-zA-Z]+', '_', sType)
            newType = '_'.join(newType.split("_"))
            if newType[-1]=="_":
                newType = newType[:-1]
            newID = dID+"_"+newType
            fout.write(sID+"\t"+newID+"\t"+dID+"\t"+sType+"\n")
        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a dictionary keyed on icgc_specimen_id, referencing its xena ID. No parms.")
    #ids()
    newIDs()
