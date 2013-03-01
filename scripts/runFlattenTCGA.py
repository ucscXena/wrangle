import os, sys, string

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

from flattenClinical import *

REALRUN =1

inDir ="data/public/TCGA/"
outDir ="data_flatten/public/TCGA/"
r = runFlatten(inDir, outDir,REALRUN,None)
