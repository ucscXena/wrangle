#!/usr/bin/env python

import argparse
import math
import re
import sys,os
import numpy
import scipy
import scipy.stats

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomicData", type=str)
    parser.add_argument("clinicalVector", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument('-p', type=str, default=None, help='optional probeMap')
    parser.add_argument('-a', type=str, default=None, help='required value of groupA in clinicalVector file')
    parser.add_argument('-b', type=str, default=None, help='required value of groupB in clinicalVector file')
    args = parser.parse_args()
    clinicalVector = readClinicalFile(args.clinicalVector)
    genomicFp = open(args.genomicData)
    if (args.p):
        probeMap = readProbeMap(args.p)
    else:
        probeMap = None
    headerTokens = genomicFp.readline().rstrip().split("\t")
    (columnsGroupA, columnsGroupB) = headerLineToGroups(headerTokens,
                                                        clinicalVector,
                                                        args.a, args.b)
    outfile = open(args.outfile, "w")
    outfile.write("%s" % (headerTokens[0]))
    if probeMap:
        outfile.write("\tgene")
    outfile.write("\tStatistic\tpValue\tMedian_%s\tMedian_%s\tDelta\n" % (args.a, args.b))

    tempOut = open(args.outfile+"_tmp", "w")
    for row in genomicFp:
        tokens = row.rstrip().split("\t")
        if len(tokens)!= len(headerTokens):
            continue
        symbol = tokens[0]
        groupATokens = [tokens[ii] for ii in columnsGroupA]
        groupBTokens = [tokens[ii] for ii in columnsGroupB]
        groupAValues = [float(xx) for xx in groupATokens if xx != "NA"]
        groupBValues = [float(xx) for xx in groupBTokens if xx != "NA"]
        (ttest, pvalue) = performTTest(groupAValues, groupBValues)
        outputLine = symbol
        if probeMap:
            if probeMap[symbol]:
                outputLine = outputLine + "\t%s" % (probeMap[symbol])
            else:
                outputLine = outputLine + "\t"
        if ttest == None or pvalue == None:
            outputLine = outputLine + "\tNA\tNA"
        else:
            outputLine= outputLine + "\t%f\t%f" % (ttest, pvalue)
        if len(groupAValues) == 0:
            outputLine = outputLine + "\tNA"
        else:
            medianA = numpy.median(groupAValues)
            outputLine = "%s\t%f" % (outputLine, medianA)
        if len(groupBValues) == 0:
            outputLine = outputLine + "\tNA"
        else:
            medianB = numpy.median(groupBValues)
            outputLine = "%s\t%f" % (outputLine, medianB)
        if len(groupAValues) == 0 or len(groupBValues) == 0:
            outputLine = outputLine + "\tNA"
        else:
            outputLine = "%s\t%f" % (outputLine, medianA - medianB)
        tempOut.write("%s\n" % (outputLine))

    tempOut.close()
    outfile.close()
    os.system("sort -nk 3 "+ args.outfile+"_tmp"+" >> "+args.outfile )
    os.system("rm "+ args.outfile+"_tmp")

def readClinicalFile(clinicalFile):
    clinicalVector = dict()
    fp = open(clinicalFile)
    fp.readline()
    for row in fp:
        tokens = row.rstrip().split("\t")
        if len(tokens) > 1:
            clinicalVector[tokens[0]] = tokens[1]
        else:
            clinicalVector[tokens[0]] = None
    fp.close()
    return(clinicalVector)

def readProbeMap(probeMapFile):
    probeMap = dict()
    fp = open(probeMapFile)
    fp.readline()
    for row in fp:
        tokens = row.rstrip().split("\t")
        probeMap[tokens[0]] = tokens[1]
    fp.close()
    return(probeMap)

def headerLineToGroups(headerTokens, clinicalVector, optionA, optionB):
    columnsGroupA = list()
    columnsGroupB = list()
    for ii in range(len(headerTokens)):
        id = headerTokens[ii]
        if clinicalVector.has_key(id):
            if clinicalVector[id] == optionA:
                columnsGroupA.append(ii)
            elif clinicalVector[id] == optionB:
                columnsGroupB.append(ii)
    return((columnsGroupA, columnsGroupB))

def performTTest(groupA, groupB):
    from scipy import special

    epsilon = 1E-06
    countA = len(groupA)
    countB = len(groupB)
    meanA = numpy.mean(groupA)
    meanB = numpy.mean(groupB)
    varA = numpy.var(groupA)
    varB = numpy.var(groupB)
    degreesOfFreedom = countA + countB - 2
    sVar = ((countA - 1) * varA + (countB - 1) * varB) / degreesOfFreedom
    if sVar == 0:
        tStat = None
        prob = None
    else:
        tStat = (meanA - meanB)/math.sqrt(sVar * (1.0/(countA + epsilon)
                                                  + 1.0/(countB + epsilon)))
        bottom = degreesOfFreedom + tStat * tStat
        prob = special.betainc(0.5*degreesOfFreedom, 0.5,
                               degreesOfFreedom/bottom)
    return(tStat, prob)




if __name__ == "__main__":
    main()
