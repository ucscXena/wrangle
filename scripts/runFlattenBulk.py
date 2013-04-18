import os, sys, string

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

from flattenClinical import *

REALRUN =0


inDir ="data/su2c/su2cPancreas/"
outDir ="data_flatten/su2c/su2cPancreas/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/su2c/su2cBreast/"
outDir ="data_flatten/su2c/su2cBreast/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/su2c/su2cEpigenetics/"
outDir ="data_flatten/su2c/su2cEpigenetics/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/gray/"
outDir ="data_flatten/gray/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/pancreas/"
outDir ="data_flatten/pancreas/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/ispy/"
outDir ="data_flatten/ispy/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/lincs/"
outDir ="data_flatten/lincs/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/public/OSHUBaylorLBL/"
outDir ="data_flatten/public/OSHUBaylorLBL/"
r = runFlatten(inDir, outDir,REALRUN,None)

inDir ="data/public/other/"
outDir ="data_flatten/public/other/"
r = runFlatten(inDir, outDir,REALRUN,None)

"""

inDir ="data/public/TCGA/"
outDir ="data_flatten/public/TCGA/"
r = runFlatten(inDir, outDir,REALRUN,None)
"""
