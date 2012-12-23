import os,sys
import optparse

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

from downloadBundle import *

PROBEMAP = "/inside/home/jzhu/cgDataJing/scripts/data_flatten/probeMap/"

parser = optparse.OptionParser()
parser.add_option("--inDir", action="store", type="string", dest="inDir")
parser.add_option("-t", action="store_true", dest="testOnly")
(options, args) = parser.parse_args()

def printUsage():
    print "python runDownloadBundle.py --inDir=inputDir -t(optional test only)\n"

    
if options.inDir==None:
    printUsage()
    sys.exit()


inDir = options.inDir
testOnly = options.testOnly

if inDir[-1]!="/":
    inDir = inDir +"/"

for dirpath, dirnames, filenames in os.walk(inDir, followlinks=True):
    if dirnames ==[]:
        print dirpath
        testVersionOnly=True
        r = downloadBundle (dirpath+"/",PROBEMAP, testVersionOnly)
        if r!=0:
            sys.exit()
if (testOnly):
    sys.exit()

for dirpath, dirnames, filenames in os.walk(inDir, followlinks=True):
    if dirnames ==[]:
        print dirpath
        testVersionOnly=False
        r = downloadBundle (dirpath+"/",PROBEMAP, testVersionOnly)
