import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def oldClin (inDir, outDir, cancer, flog,REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name
    
    #real_cancer = string.replace(cancer,"TCGA.","")
    #real_cancer = string.replace(real_cancer,".SAMPLEMAP","")
    real_cancer =cancer
    preClinMatrix = previousClin(inDir,real_cancer)
    
    preCOLs = preClinMatrix.getCOLs()
    currentCOLs = currentClin(outDir+real_cancer, real_cancer)

    
    preROWs = preClinMatrix.getROWs()
    oldClinMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_oldClin")
    oldClinMatrix.addNewRows(preROWs, {})


    for col in preCOLs:
        if col not in currentCOLs and col not in ["_PATIENT","_INTEGRATION"]:
            oldClinMatrix.addOneColWithSameValue(col,"")
            for row in preROWs:
                oldClinMatrix.setDATA(row,col,preClinMatrix.getDATA(row,col))

    if oldClinMatrix.getCOLnum() >0:
        output (outDir+real_cancer+"/", oldClinMatrix, real_cancer)


def output (dir, oldClinMatrix, cancer):
    #data output
    name = oldClinMatrix.getName()
    fout=open(dir+name,'w')
    oldClinMatrix.store(fout,True)
    fout.close()

    #JSON
    J={}
    J["cgDataVersion"]=1
    J["redistribution"]= True
    J["dataProducer"]= "UCSC"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"]= "old clinical data that is deprecated in the current version"
    J["name"]= name
    J["type"]= "clinicalMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #cFjson
    cFJ={}
    cFJ["name"]=J["name"]+"_clinFeat"
    cFJ["type"]="clinicalFeature"

    J[":clinicalFeature"] = cFJ["name"]

    #data JSON output
    fout=open(dir+name+".json",'w')
    fout.write( json.dumps( J, indent=-1 ) )
    fout.close()

    #clinicalFeature
    cFfile =dir+ cFJ["name"]
    fout = open(cFfile,"w")
    cols = oldClinMatrix.getCOLs()
    for col in cols:
        fout.write( col+"\tshortTitle\t_DEPRECATED_FEATURE: "+col+"\n")
        fout.write( col+"\tlongTitle\t_DEPRECATED_FEATURE: "+col+"\n")
    fout.close()
    
    #clinicalFeature JSON output
    fout=open(dir+cFJ["name"]+".json",'w')
    fout.write( json.dumps( cFJ, indent=-1 ) )
    fout.close() 
    
def previousClin (dir, cancer):
    file = dir +cancer+"_clinicalMatrix"
    preClinMatrix = ClinicalMatrixNew(file, "previous")
    return preClinMatrix

def currentClin (dir,cancer):
    ignore=0
    bookDic=cgWalk(dir,ignore)
    sampleMaps = collectSampleMaps(bookDic)
    allMaps = sampleMaps.keys()
    if len(allMaps)!=1:
        print "ERROR"
        return 0
    sampleMap=allMaps[0]

    currentFeatures=[]
    datasets = collectNamesBelongToSampleMap(bookDic, sampleMap)
    for name in datasets:
        obj= bookDic[name]
        if obj['type']!="clinicalMatrix":
            continue
        if name == "clinical_"+cancer+"_oldClin":
            continue

        path = obj['path']
        fin = open(path,'r')
        features = string.split(string.strip(fin.readline()),"\t")[1:]
        for feature in features:
            if feature not in currentFeatures:
                currentFeatures.append(feature)

    return currentFeatures


    
