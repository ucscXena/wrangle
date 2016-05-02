import os, sys

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *

def autoStateOrder(clinMatrix, clinFeature):
    features = clinMatrix.getCOLs()

    for feature in features:
        if clinMatrix.isTypeCategory(feature):
            states = clinMatrix.getColStates(feature)
            if len(states)==1:
                continue
            if len(states)>100:
                continue
            if clinFeature.getStateOrder(feature):
                continue
            states.sort()
            stateOrder = states
            clinFeature.setFeatureValueType(feature, "category")
            clinFeature.setFeatureStates(feature, states)
            clinFeature.setFeatureStateOrder(feature, stateOrder)

if len(sys.argv[:])!=4:
    print "python autoStateOrderFeatureFile.py clinMatrix clinFeatureOut clinFeatureIn"
    sys.exit()

clinMatrixFile = sys.argv[1]
clinFeatureFile = sys.argv[3]
output = sys.argv[2]

clinMatrix =ClinicalMatrixNew(clinMatrixFile,"matrix")

if os.path.exists(clinFeatureFile):
    clinFeature = ClinicalFeatureNew(clinFeatureFile,'feature')
    autoStateOrder(clinMatrix, clinFeature)
    fout = open(output,'w')
    clinFeature.store(fout)
else:
    clinFeature = None

