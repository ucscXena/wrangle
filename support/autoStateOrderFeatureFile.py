import os, sys, json

sys.path.insert(0, os.path.dirname(sys.argv[0])+"/../CGDataNew")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *

def autoStateOrder(clinMatrix, clinFeature):
    features = clinMatrix.getCOLs()

    for feature in features:
        if clinMatrix.isTypeCategory(feature):
            states = clinMatrix.getColStates(feature)
            if len(states)<=1:
                continue
            if len(states)>100:
                continue
            if clinFeature.getStateOrder(feature):
                continue
            print (feature)
            states.sort()
            stateOrder = states
            clinFeature.setFeatureValueType(feature, "category")
            clinFeature.setFeatureStates(feature, states)
            clinFeature.setFeatureStateOrder(feature, stateOrder)

if len(sys.argv[:]) <3:
    print ("python autoStateOrderFeatureFile.py clinMatrix clinFeatureOut clinFeatureIn(optional)")
    print 
    sys.exit(1)

clinMatrixFile = sys.argv[1]
output = sys.argv[2]

clinMatrix =ClinicalMatrixNew(clinMatrixFile, "matrix")
clinFeature = None
print (clinMatrix.getCOLs())
if len(sys.argv[:])==4:
    clinFeatureFile = sys.argv[3]
    if os.path.exists(clinFeatureFile):
        clinFeature = ClinicalFeatureNew(clinFeatureFile,'feature')
    else:
        print (sys.argv[3],"does not exist")
        sys.exit(1)
else:
    clinFeature = ClinicalFeatureNew(None,'feature')
    fout = open(output+".json",'w')
    J={}
    J["type"]="clinicalFeature"
    fout.write(json.dumps(J, indent=2))
    fout.close()
    
autoStateOrder(clinMatrix, clinFeature)
fout = open(output,'w')
clinFeature.store(fout)
