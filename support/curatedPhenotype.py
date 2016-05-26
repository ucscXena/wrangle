import os, sys,json,getopt

sys.path.insert(0, os.path.dirname(sys.argv[0])+"/../CGDataNew")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *

def getCuratedPhenotype():
    dic ={
        "_primary_disease":{
            "shortTitle":"primary_disease",
            "longTitle":"primary_disease",
            "type":"category"
        },
        "_primary_site":{
            "shortTitle":"primary_site",
            "longTitle":"primary_site",
            "type":"category"
        },
        "_cohort":{
            "shortTitle":"cohort",
            "longTitle":"cohort",
            "type":"category"
        },
        "_gender":{
            "shortTitle":"gender",
            "longTitle":"gender",
            "type":"category"
        },
        "_sample_type":{
            "shortTitle":"sample type",
            "longTitle":"sample type",
            "type":"category"
        },
        "_OS_IND":{
            "shortTitle":"overall survival indicator",
            "longTitle":"overall survival indicator",
            "type":"category"
        },
        "_OS_UNIT":{
            "shortTitle":"overall survival time unit",
            "longTitle":"overall survival time unit",
            "type":"category"
        },
        "_OS":{
            "shortTitle":"overall survival time",
            "longTitle":"overall survival time",
            "type":"float"
        },
        "_age_at_diagnosis":{
            "shortTitle":"age_at_diagnosis",
            "longTitle":"age_at_diagnosis",
        }
    }
    return dic

def curatedPhenotypeClinFeature (clinFeature):
    curated = getCuratedPhenotype()
    for feature in curated.keys():
        data = curated[feature]
        shortTitle = data["shortTitle"]
        longTitle = data["longTitle"]
        clinFeature.setFeatureShortTitle(feature, shortTitle)
        clinFeature.setFeatureLongTitle(feature, longTitle)
        if "type" in data:
            featureType = data["type"]
            clinFeature.setFeatureValueType(feature, featureType)
    return


#http://www.tutorialspoint.com/python/python_command_line_arguments.htm
#https://docs.python.org/3.1/library/getopt.html
try:
    opts, args = getopt.getopt(sys.argv[1:],"",["run"])
except getopt.GetoptError:
    print "python curatedPhenotype.py originalClinFeature(optional) --run"
    sys.exit()


output = "newClinFeature"
clinFeature = None

if len(args)!=0:
    clinFeatureFile = args[0]
    if os.path.exists(clinFeatureFile):
        clinFeature = ClinicalFeatureNew(clinFeatureFile,'feature')
    else:
        print args[0],"does not exist"
        sys.exit()
else:
    clinFeature = ClinicalFeatureNew(None,'feature')
    fout = open(output+".json",'w')
    J={}
    J["type"]="clinicalFeature"
    fout.write(json.dumps(J, indent=2))
    fout.close()

curatedPhenotypeClinFeature(clinFeature)

fout = open(output,'w')
clinFeature.store(fout)
print "output:", output
