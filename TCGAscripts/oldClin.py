import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def oldClin (inDir, outDir, cancer, flog,REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    currentClinMatrix = currentClin(outDir+cancer, cancer)
    preClinMatrix = previousClin(inDir,cancer)

    preCOLs = preClinMatrix.getCOLs()
    currentCOLs = currentClinMatrix.getCOLs()

    preROWs = preClinMatrix.getROWs()
    oldClinMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_oldClin")
    oldClinMatrix.addNewRows(preROWs, {})

    for col in preCOLs:
        if col not in currentCOLs:
            oldClinMatrix.addOneColWithSameValue(col,"")
            for row in preROWs:
                oldClinMatrix.setDATA(row,col,preClinMatrix.getDATA(row,col))

    if oldClinMatrix.getCOLnum() >0:
        output (outDir+cancer+"/", oldClinMatrix, cancer)

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
    file =dir+cancer+"_clinicalMatrix"
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
    finalClinMatrix= ClinicalMatrixNew(None,"tmp")
    datasets = collectNamesBelongToSampleMap(bookDic, sampleMap)
    datasetsOrdered =[] #only the ClinicalMatrix ordered list
    for name in datasets:  
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            if obj.has_key('outOfDate') and obj['outOfDate'] in ["yes", "Yes","YES"]:
                datasetsOrdered.append(name)
            elif not obj.has_key('outOfDate') and  not obj.has_key('upToDate'):
                datasetsOrdered.insert(0,name)

    upToDateSets={}
    for name in datasets:  
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            if obj.has_key('upToDate') :
                upToDateSets[obj['upToDate']]=name
    keys= upToDateSets.keys()
    keys.sort()
    for version in keys:
        name = upToDateSets [version]
        datasetsOrdered.insert(0,name)

    for name in datasetsOrdered:
        if name == "clinical_"+cancer+"_oldClin":
            continue
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            #get matrix obj
            path = obj['path']
            name = obj['name']
            cMatrix = ClinicalMatrixNew(path,name)
            
            if finalClinMatrix==None:
                finalClinMatrix= cMatrix
                
            #merge final and cMatrix
            if finalClinMatrix != cMatrix:
                r = finalClinMatrix.addNewCols(cMatrix,validation=True)
                if r!=True:
                    print "Fail to merge"
                    return False
    return finalClinMatrix


    
