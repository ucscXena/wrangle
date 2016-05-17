import sys,string, os
import config

thisDir =  os.path.dirname(sys.argv[0]);
if thisDir =="":
    thisDir ="."

sys.path.insert(0, thisDir +"/../CGDataNew")
from ClinicalMatrixNew import *

def modifyClinMatrix(clinMatrix, col, keyCol, dataDic):
    samples = clinMatrix.getROWs()
    for sample in samples:
        key = clinMatrix.getDATA(sample,keyCol)
        if (key  and  (key in dataDic)):
            value = dataDic[key]
            clinMatrix.setDATA(sample, col, value)

if len(sys.argv[:])!=3:
    print "python addPhenotypeByAPI.py inputClinicalMatrix outputClinicalMatrix"
    sys.exit()

config.getPrimaryDisease()

inFile = sys.argv[1]
outFile = sys.argv[2]
clinMatrix =ClinicalMatrixNew(inFile,"matrix")

key_col = "project_code"
feature = "_primary_site"
dataDic = config.getPrimarySite()
clinMatrix.addOneColWithSameValue(feature,"")
modifyClinMatrix(clinMatrix, feature, key_col, dataDic)

feature = "_primary_disease"
dataDic = config.getPrimaryDisease()
clinMatrix.addOneColWithSameValue(feature,"")
modifyClinMatrix(clinMatrix, feature, key_col, dataDic)

fout = open(outFile,'w')
clinMatrix.store(fout)
fout.close()
