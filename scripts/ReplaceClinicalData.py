import string, os, sys,stat
import json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

from CGDataLib import *

if len(sys.argv[:]) != 4:
    print "python ReplaceClinicalData.py clinFileIn codeJson clinFileout"
    print
    sys.exit()

oldMatrix = ClinicalMatrixNew(sys.argv[1], "old")

fin = open(sys.argv[2],'r')
J = json.loads(fin.read())
fin.close()

for feature in oldMatrix.getCOLs():
    if J.has_key(feature):
        print feature
        if J[feature].has_key("code"):
            code= J[feature]["code"]
            for key in code:
                old =key
                new= code[key]
                oldMatrix.replaceValueInCol(feature, old, new)

fout=open(sys.argv[3],'w')
oldMatrix.store(fout,True)
fout.close()

        
