import string, os, sys
import json,datetime

PATHPATTERN="oneoff"

import TCGAUtil
sys.path.insert(0,"../CGDataNew/")
from ClinicalMatrixNew import *
from ClinicalFeatureNew import *
from CGDataUtil import *

#/data/TCGA/tcgaDataOneOff/*/

def OneOffClinical(inDir, outDir, cancer,flog,REALRUN):
    if cancer in ["COAD","READ"]:
        deriveCancer="COADREAD"
        process (inDir,outDir,deriveCancer,flog,PATHPATTERN,cancer)
    process (inDir,outDir,cancer,flog,PATHPATTERN,cancer)

def process (inDir,outDir,cancer,flog,PATHPATTERN,originCancer):
    #print status
    print cancer, __name__

    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    for file in os.listdir(inDir):
        clinMatrix = None
        clinFeature =None
        clinFfile=""

        #find the file
        #clinMatrix

        if file[0:6]== PATHPATTERN and os.path.exists(inDir+ file+".json") :
            pass
        else:
            continue

        infile = inDir+file

        #json file processing (validation)
        fjson= open(infile+".json","U")
        J =json.load(fjson)
        fjson.close()

        if J["type"]!="clinicalMatrix":
            continue

        #clinFeature
        if J.has_key(":clinicalFeature"):
            clinFname = J[":clinicalFeature"]
       
            for clinFfile in os.listdir(inDir):
                #find the file
                if not os.path.exists(inDir+ clinFfile+".json"):
                    continue

                fjson= open(inDir+clinFfile+".json","U")
                clinFJ =json.load(fjson)
                fjson.close()

                #data processing
                if clinFJ["type"]=="clinicalFeature" and clinFJ["name"]==clinFname:
                    print originCancer, cancer
                    if cancer != originCancer:
                        clinFname= clinFname+"_"+originCancer
                        clinFJ["name"]=clinFname
                    clinFeature= ClinicalFeatureNew(inDir+clinFfile,clinFname)
                    for feature in clinFeature.getFeatures():
                        if TCGAUtil.featurePriority.has_key(cancer):
                            if TCGAUtil.featurePriority[cancer].has_key(feature):
                                priority= TCGAUtil.featurePriority[cancer][feature]
                                clinFeature.setFeaturePriority(feature, priority)
                                clinFeature.setFeatureVisibility(feature, "on")
                    break
                
        #data processing
        clinMatrix = ClinicalMatrixNew(infile, J["name"], False, clinFeature)
        clinMatrix.removeCols(["ethnicity","race","jewish_origin"])
        clinMatrix.replaceValue("null","")
        clinMatrix.replaceValue("NULL","")
        clinMatrix.replaceValue("Null","")
        clinMatrix.replaceValue("NA","")
        clinMatrix.replaceValue("[null]","")
        clinMatrix.replaceValue("[NULL]","")
        clinMatrix.replaceValue("[Null]","")
        clinMatrix.replaceValue("[NA]","")
        clinMatrix.replaceValue("[Not Available]","")
        clinMatrix.replaceValue("[Not Reported]","")
        clinMatrix.replaceValue("[Not Applicable]","")
        clinMatrix.replaceValue("[Not Requested]","")
        clinMatrix.replaceValue("[Completed]","")
        clinMatrix.replaceValue("[Pending]","")
        clinMatrix.replaceValue("[]","")
        clinMatrix.replaceValue("Not Tested","")

        if cancer != originCancer:
            clinMatrix.addOneColWithSameValue("cancer type",originCancer)

        #json file processing (validation)
        fjson= open(infile+".json","U")
        J =json.load(fjson)
        fjson.close()
        if cancer != originCancer:
            J['name'] = J['name'] +"_"+originCancer
            J[":sampleMap"]="TCGA."+cancer+".sampleMap"
            
        name = trackName_fix(J['name'])
        if name ==False:
            message = "bad object name, need fix otherwise break loader, too long "+J["name"]
            print message
            flog.write(message+"\n")
            return
        else:
            J["name"]=name

        if cancer != originCancer and J.has_key(":clinicalFeature"):
            J[":clinicalFeature"] =  J[":clinicalFeature"] +"_"+originCancer

        J["cgDataVersion"]=1

        #output matrix
        if cancer != originCancer:
            outfile = outDir+cancer+"/"+file+"_"+originCancer
        else:
            outfile = outDir+cancer+"/"+file
            
        oHandle = open(outfile,"w")
        clinMatrix.store(oHandle, validation=True)
        oHandle.close()

        fjson = open(outfile+".json","w")
        json.dump(J, fjson, indent=-1)
        fjson.close()

        #output clinFeature 
        if clinFeature:
            if cancer != originCancer:
                outfile = outDir+cancer+"/"+clinFfile+"_"+originCancer
            else:
                outfile = outDir+cancer+"/"+clinFfile
            fout=open(outfile,'w')
            clinFeature.store(fout)
            fout.close()

            clinFJ["cgDataVersion"]=1
            fjson = open(outfile+".json","w")
            json.dump(clinFJ, fjson, indent=-1)
            fjson.close()
    return

