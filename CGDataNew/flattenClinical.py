import sys,os,string,copy,json
from SampleMapNew import *
from CGDataLib import *
from ClinicalMatrixNew  import *
from ClinicalFeatureNew  import *
from IntegrationId  import *
from CGDataUtil import *

def runFlatten(inDir, outDir,REALRUN, onlyGenomicSamples, SMAPNAME=None):
    dir = inDir
    bookDic={}
    sampleMaps={}
    ignore=0
    bookDic=cgWalk(dir,ignore)
    if not bookDic :
        print "repo has problem"
        return 0
    sampleMaps = collectSampleMaps(bookDic)
    missingMaps= collectMissingSampleMaps(bookDic)

    allMaps = sampleMaps.keys()
    allMaps.extend(missingMaps.keys())

    for sampleMap in allMaps:
        if SMAPNAME and SMAPNAME!=sampleMap:
            print "skip", sampleMap
            continue

        print sampleMap
        path = bookDic[sampleMap]['path']
        if os.path.abspath(path) in [ \
            "/inside/home/jzhu/cgDataJing/scripts/data/public/TCGA/PANCAN/TCGA.PANCAN.sampleMap", \
                "/inside/home/jzhu/cgDataJing/scripts/data/public/TCGA/PANCAN12/TCGA.PANCAN12.sampleMap" ]:
            print "ignore "+path
            continue

        if sampleMap in missingMaps:
            #construct an empty sampleMap
            sMap = SampleMapNew(None,sampleMap)
            #fill sMap with individual nodes, no connection
            changed = checkIdsAllIn (sMap, bookDic)
            #build connection
        else:
            name = bookDic[sampleMap]['name']
            fin = open(path,'r')
            sMap = SampleMapNew(fin,name)
            if not sMap.getName():
                print "Fail to initiate", name
                return 0
            fin.close()
            changed = checkIdsAllIn (sMap, bookDic)
        
        if REALRUN in [0,1]:
            r = flattenEachSampleMap(sMap,bookDic,onlyGenomicSamples)
            if r== False:
                return 0
            finalClinicalMatrix,finalClinicalMatrixJSON,finalClinFeature,finalClinFeatureJSON= r
            if finalClinicalMatrix.getROWnum()!=0:
                outputEachSampleMapRelated(outDir, bookDic, sMap,
                                           finalClinicalMatrix,finalClinicalMatrixJSON,
                                           finalClinFeature,finalClinFeatureJSON,REALRUN)
        if REALRUN == -2:
            finalClinFeature = flattenForClinicalFeature(sMap, bookDic)
            outputForClinFeature(outDir,sMap, finalClinFeature)
            
        cpGenomicEachSample(REALRUN, outDir, bookDic, sMap)
        cpProbeMaps(REALRUN,outDir,bookDic,sMap)
            
    return 1

def flattenForClinicalFeature(sMap, bookDic):
    clinFeatures=[]
    finalClinFeature=None
    sampleMap = sMap.getName()
    datasets = collectNamesBelongToSampleMap(bookDic, sampleMap)
    for name in datasets:  
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            clinFeature=None
            #clinFeature obj
            if obj.has_key(':clinicalFeature'):
                path=  bookDic[obj[':clinicalFeature']]['path']
                neme = bookDic[obj[':clinicalFeature']]['name']
                clinFeature = ClinicalFeatureNew(path, name)

            #get matrix obj
            path = obj['path']
            name = obj['name']
            cMatrix = ClinicalMatrixNew(path,name,False, clinFeature)

            if clinFeature:
                clinFeatures.append(clinFeature)

    fout = open(".tmp",'w')
    fout.close()
    for clinF in clinFeatures:
        fout = open(".tmptmp",'w')
        clinF.store(fout)
        fout.close()
        os.system("cat .tmptmp >> .tmp")
    fin = open(".tmp",'r')
    jsonName=  trackName_fix(sampleMapBaseName(sMap)+"_clinicalFeature")
    finalClinFeature =ClinicalFeatureNew(fin,jsonName)
    if not finalClinFeature.isValid():
        print "final clinFeature file .tmp is invalid"
        return 0
    fin.close()

    #vis exceptions
    VIS_limit=4
    if bookDic.has_key(sMap.getName()) and bookDic[sMap.getName()].has_key("VIS"):
        VIS_limit= bookDic[sMap.getName()]["VIS"]
    finalClinFeature.fillInPriorityVisibility(VIS_limit)
    finalClinFeature.setFeatureShortTitle("_PATIENT","_PATIENT_ID")
    finalClinFeature.setFeatureLongTitle("_PATIENT","_PATIENT_ID")
    finalClinFeature.setFeatureValueType("_PATIENT","category")
    finalClinFeature.setFeatureShortTitle("_INTEGRATION","_SAMPLE_ID")
    finalClinFeature.setFeatureLongTitle("_INTEGRATION","_SAMPLE_ID")
    finalClinFeature.setFeatureValueType("_INTEGRATION","category")
    return finalClinFeature


def flattenEachSampleMap(sMap, bookDic,onlyGenomicSamples):
    sampleMap = sMap.getName()

    jsonName= trackName_fix(sampleMapBaseName(sMap)+"_clinicalMatrix")
    finalClinMatrix= ClinicalMatrixNew(None,jsonName)
    finalClinMatrixJSON={}
    finalClinMatrixJSON["name"]=jsonName
    finalClinMatrixJSON["type"]="clinicalMatrix"
    finalClinMatrixJSON["path"]=""
    finalClinMatrixJSON[":sampleMap"]=sampleMap

    clinFeatures=[]
    finalClinFeatureJSON=None
    finalClinFeature=None

    # add all ids to sMap
    sMapChanged= checkIdsAllIn(sMap, bookDic)

    #build initial clinical Matrix with sampleMap ids, all with empty data
    emptyData={}
    success = finalClinMatrix.addNewRows(sMap.getNodes(),emptyData)
    if not success:
        print "fail to add all initial ids from sampleMap"

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
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            clinFeature=None
            #clinFeature obj
            if obj.has_key(':clinicalFeature'):
                path=  bookDic[obj[':clinicalFeature']]['path']
                neme = bookDic[obj[':clinicalFeature']]['name']
                clinFeature = ClinicalFeatureNew(path, name)

            #get matrix obj
            path = obj['path']
            name = obj['name']

            cMatrix = ClinicalMatrixNew(path,name,False, clinFeature)
            
            if finalClinMatrix==None:
                finalClinMatrix= cMatrix
                
            if finalClinMatrixJSON==None:
                finalClinMatrixJSON= obj

            #merge final and cMatrix
            if finalClinMatrix != cMatrix:
                print "name=",cMatrix.getName()
                r = finalClinMatrix.addNewCols(cMatrix,validation=True)
                if r!=True:
                    print "Fail to merge"
                    return False

            #add clinFeature
            if clinFeature:
                clinFeatures.append(clinFeature)
            
            #merge finalClinMatrixJSON with new json
            if finalClinMatrixJSON != obj:
                jsonName= trackName_fix(sampleMapBaseName(sMap)+"_clinicalMatrix")
                finalClinMatrixJSON= cgDataMergeJSON(finalClinMatrixJSON, obj, jsonName)

            # final ClinFeature json
            if clinFeature:
                clinFeatureJSON = bookDic[obj[':clinicalFeature']]
                if finalClinFeatureJSON==None:
                    finalClinFeatureJSON= clinFeatureJSON
                else:
                    jsonName=  trackName_fix(sampleMapBaseName(sMap)+"_clinicalFeature")
                    finalClinFeatureJSON["version"]=datetime.date.today().isoformat() 
                    finalClinFeatureJSON["type"]="clinicalFeature"
                    finalClinFeatureJSON["name"]=jsonName
                    finalClinFeatureJSON["cgDataVersion"]=1

    #final clinicalFeature
    if finalClinFeatureJSON:
        fout = open(".tmp",'w')
        fout.close()
        for clinF in clinFeatures:
            fout = open(".tmptmp",'w')
            clinF.store(fout)
            fout.close()
            os.system("cat .tmptmp >> .tmp")
        fin = open(".tmp",'r')
        finalClinFeature =ClinicalFeatureNew(fin,finalClinFeatureJSON['name'])
        if not finalClinFeature.isValid():
            print "final clinFeature file .tmp is invalid"
            return 0
        fin.close()
    
    #SURVIVAL analysis data
    foundE=0
    foundT=0
    if finalClinFeature:
        features= finalClinFeature.getFeatures()
        for feature in features:
            sameAs = finalClinFeature.getFeatureSameAs(feature)
            if sameAs =="_TIME_TO_EVENT":
                #check there is only one parameter is set to be _TIME_TO_EVENT
                if foundT==1:
                    print "ERROR there is already _TIME_TO_EVENT"
                    continue
                #check matrix does not have _TIME_TO_EVNET
                if sameAs in finalClinMatrix.getCOLs():
                    print "ERROR there is already _TIME_TO_EVENT in matrix"
                    continue
                #data check need to check these are floats or "" in both clinFeature and clinMatrix 
                GOOD=1
                if finalClinMatrix.isTypeFloat(feature)!= True:
                    print "ERROR _TIME_TO_EVENT parent feature values are not correct", finalClinMatrix.getColStates(feature)
                    GOOD=0
                if GOOD:
                    foundT=1
                    finalClinMatrix.addNewColfromOld(sameAs, feature)
                    finalClinFeature.setFeatureValueType(sameAs,"float")
                    
            if sameAs =="_EVENT":
                #check there is only one parameter is set to be _EVENT
                if foundE==1:
                    print "ERROR there is already _EVENT"
                    continue
                #check matrix does not have _EVNET
                if sameAs in finalClinMatrix.getCOLs():
                    print "ERROR there is already _EVENT in matrix"
                    continue
                #data check
                GOOD=1
                states= finalClinMatrix.getColStates(feature)
                """
                for state in states:
                    if state not in [0,1,"0","1",""]:
                        print "ERROR _EVENT values are not correct", state
                        GOOD=0
                        break
                """
                if len(states) not in [2,3]:
                    GOOD=0
                if len(states)==3 and states.count('')!=1:
                    GOOD=0

                if GOOD:
                    foundE=1
                    finalClinMatrix.addNewColfromOld(sameAs, feature)
                    finalClinFeature.setFeatureValueType(sameAs,"category")
                    #finalClinFeature.setFeatureStates(sameAs,["0","1"])
                    #finalClinFeature.setFeatureStateOrder(sameAs,["0","1"])

    #clinical data push down
    roots = sMap.allRoots()
    for root in roots:
        r = finalClinMatrix.pushToChildren (root,sMap)
        if r != True:
            print "Fail to push down"
            return 0
    print "after clinical push down", sampleMap,finalClinMatrix.getROWnum()
    
    # collect all genomic data
    keepSamples  = getAllGenomicIds(sMap, bookDic)

    # removing rows without genomic data from  clinical data matrix due to mysql enum limitation
    # should remove this step after the display functionality is done better, currently cgb clinical data range in feature control panel shows the full range of clinical data without checking if the specific track/dataset has the full value range.
    if onlyGenomicSamples:
        print "genomic sample count", len(keepSamples)
        success= finalClinMatrix.onlyKeepRows(keepSamples)
        if not success:
            print "fail to remove extra rows"
        else:
            print "after keeping sample with genomic data", finalClinMatrix.getROWnum()
        
    #add to the clinical matrix any samples with genomic data but no clinical data
    emptyData={}
    for col in finalClinMatrix.getCOLs():
        emptyData[col]=""
    success = finalClinMatrix.addNewRows(keepSamples,emptyData)
    if not success:
        print "fail to add new roows"
    else:
        print "after adding all genomic data", finalClinMatrix.getROWnum()

    if finalClinMatrix.validate() != True:
        print "Fail to validate"
        cMatrix = oldCMatrix
        return 0
    # end of collecting all genomic data
    
    #code to remove blacklist samples and all its descendants

    badList= badListSelfAndDescendants (sMap, bookDic)
    if badList!=[]:
        #remove badList samples
        finalClinMatrix.removeRows(badList, True)
        print "after remove badList", finalClinMatrix.getROWnum()
        
    #identify empty features 
    badFeatures= finalClinMatrix.findBadColsNotRemove()
    finalBadFeatures=[]

    if finalClinFeature:
        for feature in badFeatures:
            #get short label
            shortTitle = finalClinFeature.getShortTitle(feature)
            if not shortTitle:
                print feature,"remove"
                finalBadFeatures.append(feature)
            elif shortTitle[:19]!="_DEPRECATED_FEATURE":
                finalBadFeatures.append(feature)
            else:
                print shortTitle,"not remove"
    else:
        finalBadFeatures =badFeatures[:]
        
    #remove bad features
    finalClinMatrix.removeCols(finalBadFeatures)
    
    if badFeatures:
        print "remove features", badFeatures

    # add _PATIENT col
    if finalClinMatrix.addColRoot(sMap) == None:
        print "Fail to addColRoot"
        return 0
            
    # add _INTEGRATION col
    intList=[]
    if bookDic.has_key(sampleMap) and bookDic[sampleMap].has_key(":integrationId"):
        intName=bookDic[sampleMap][":integrationId"]
        fin= open(bookDic[intName]["path"],"r")
        intId = IntegrationId(intName,fin)
        intList = intId.getList()
    finalClinMatrix.addColIntegration(sMap,intList)
    
                
    # final ClinFeature json
    if finalClinFeatureJSON==None:
        jsonName=  trackName_fix(sampleMapBaseName(sMap)+"_clinicalFeature")
        finalClinFeatureJSON= {}
        finalClinFeatureJSON["version"]=datetime.date.today().isoformat() 
        finalClinFeatureJSON["type"]="clinicalFeature"
        finalClinFeatureJSON["cgDataVersion"]=1
        finalClinFeatureJSON["name"]=jsonName
        finalClinFeatureJSON["path"]=""
        finalClinFeature = ClinicalFeatureNew (None, finalClinFeatureJSON["name"])

    #final clinicalFeature
    if finalClinFeature:
        finalClinFeature.removeFeatures(finalBadFeatures)
        finalClinFeature.cleanState()
        finalClinFeature.checkFeatureWithMatrix(finalClinMatrix)
        #clinicalFeature fillin ValueType
        finalClinFeature.fillInValueTypeWithMatrix(finalClinMatrix)
        #clinicalFeature fillin missing features
        finalClinFeature.fillInFeaturesWithMatrix(finalClinMatrix)
        #clinicalFeature fillin short and long titles
        finalClinFeature.fillInTitles()
        #clinicalFeature fillin priority visibility

        #vis exceptions
        VIS_limit=4
        if bookDic.has_key(sMap.getName()) and bookDic[sMap.getName()].has_key("VIS"):
            VIS_limit= bookDic[sMap.getName()]["VIS"]
        finalClinFeature.fillInPriorityVisibility(VIS_limit)
        
        finalClinFeature.setFeatureShortTitle("_PATIENT","_PATIENT_ID")
        finalClinFeature.setFeatureLongTitle("_PATIENT","_PATIENT_ID")
        finalClinFeature.setFeatureValueType("_PATIENT","category")
        finalClinFeature.setFeatureShortTitle("_INTEGRATION","_SAMPLE_ID")
        finalClinFeature.setFeatureLongTitle("_INTEGRATION","_SAMPLE_ID")
        finalClinFeature.setFeatureValueType("_INTEGRATION","category")
        
    print sampleMap,finalClinMatrix.getROWnum()
    return finalClinMatrix,finalClinMatrixJSON, finalClinFeature, finalClinFeatureJSON

def outputForClinFeature(outDir,sMap, finalClinFeature):
    dataPackageDir = outDir + sampleMapBaseName(sMap)
    #output clinicalFeature data by concatenation
    if not os.path.exists( dataPackageDir ):
        os.makedirs( dataPackageDir )
    targetfile = dataPackageDir+"/"+sampleMapBaseName(sMap)+"_clinicalFeature"
    fout = open(targetfile,'w')
    finalClinFeature.store(fout)
    fout.close()
        
def outputEachSampleMapRelated(outDir, bookDic, sMap,
                               finalClinMatrix,finalClinMatrixJSON,
                               finalClinFeature,finalClinFeatureJSON,REALRUN):
    if REALRUN not in [0,1]:
        return
    
    sampleMap = sMap.getName()
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    dataPackageDir = outDir + sampleMapBaseName(sMap)
    if not os.path.exists( dataPackageDir ):
        os.makedirs( dataPackageDir )
    else:  #clean out the existing dataPackageDir to remove all old datafiles if realrun
        if REALRUN:
            os.system("rm -rf "+dataPackageDir+"/*")
        
    #copy sampleMap and .json
    path = bookDic[sampleMap]['path']
    os.system("cp "+path+".json "+dataPackageDir+"/")
    fout = open(dataPackageDir+"/"+os.path.basename(path),'w')
    sMap.store(fout)
    fout.close()

    #copy integratiionId and .json, ignore it for now
    """
    if bookDic[sampleMap].has_key(":integrationId"):
        intName=bookDic[sampleMap][":integrationId"]
        path = bookDic[intName]['path']
        os.system("cp "+path+" "+dataPackageDir+"/")
        os.system("cp "+path+".json "+dataPackageDir+"/")
    """
        
    #output clinicalMatrix data
    fout = open(dataPackageDir+"/"+sampleMapBaseName(sMap)+"_clinicalMatrix",'w')
    finalClinMatrix.store(fout)
    fout.close()
        
    #copy clincalMatrix json
    finalClinMatrixJSON.pop('path')
    fout = open(dataPackageDir+"/"+sampleMapBaseName(sMap)+"_clinicalMatrix.json",'w')
    if finalClinFeatureJSON != None:
        finalClinMatrixJSON[':clinicalFeature']= finalClinFeatureJSON["name"]
    if not finalClinMatrixJSON.has_key("version"):
        finalClinMatrixJSON["version"]= bookDic[sampleMap]["version"]
    fout.write(json.dumps(finalClinMatrixJSON, indent=-1))
    fout.close()

    #output clinicalFeature data by concatenation
    if finalClinFeatureJSON != None:
        targetfile = dataPackageDir+"/"+sampleMapBaseName(sMap)+"_clinicalFeature"
        fout = open(targetfile,'w')
        finalClinFeature.store(fout)

    #output clinicalFeatureJSON
    if finalClinFeatureJSON != None:
        finalClinFeatureJSON.pop('path')
        fout = open(dataPackageDir+"/"+sampleMapBaseName(sMap)+"_clinicalFeature.json",'w')
        fout.write(json.dumps(finalClinFeatureJSON, indent=-1))
        fout.close()

def cpGenomicEachSample(REALRUN, outDir, bookDic, sMap):
    if REALRUN==-2:
        return
    
    dataPackageDir = outDir + sampleMapBaseName(sMap)
    for name in  collectNamesBelongToSampleMap(bookDic, sMap.getName()):
        J = bookDic[name]
        if J['type'] not in [ "clinicalMatrix","sampleMap"]:
            path = J['path']
            if not J.has_key("redistribution"):
                J["redistribution"]=True
            elif J["redistribution"] in ["false","FALSE","False", False,"0", 0 ,"no","NO","No"]:
                J["redistribution"]=False
            elif J["redistribution"] in ["true","TRUE","True", True,"1", 1 ,"yes","YES","Yes"]:
                J["redistribution"]=True

            # check security tag, if it is private (any forms) set to "private", if not present, set to "public"
            if not J.has_key("security"):
                J["security"]="public"
            elif J["security"] in ["private","PRIVATE","Private","controlled-access"]:
                J["security"]="private"
                
            outfile = dataPackageDir+"/"+os.path.basename(path+".json")
            oHandle=open(outfile,"w")
            oHandle.write( json.dumps( J, indent=-1 ) )
            oHandle.close()

            if REALRUN == 1:
                #code to remove blacklist samples and all its descendants
                badList= badListSelfAndDescendants (sMap, bookDic)
                if badList==[]:
                    os.system("cp "+path+" "+dataPackageDir+"/")
                else:
                    if J['type']=="genomicMatrix":
                        goodCols=["1"]
                        badCols=[]
                        fin =open (path,'r')
                        samples = string.split(string.strip(fin.readline()),"\t")
                        fin.close()
                        for i in range(1,len(samples)):
                            if samples[i] not in badList:
                                goodCols.append(str(i+1))
                            else:
                                badCols.append(str(i+1))
                        if badCols==[]:
                            os.system("cp "+path+" "+dataPackageDir+"/")
                        else:
                            #print badCols
                            s = string.join(goodCols,",")
                            os.system("cut -f "+s+" "+path+" > "+dataPackageDir+"/"+os.path.basename(path))
                    if J['type'] in ["genomicSegment","mutationVector"]:
                        s = "grep -v "+badList[0]+" "+path
                        for i in range (1,len(badList)):
                            s= s+" | grep -v "+ badList[i]
                        s =s +" > " +dataPackageDir+"/"+os.path.basename(path)
                        os.system(s)
    return

def badListSelfAndDescendants (sMap, bookDic):
    #code to remove blacklist samples and all its descendants
    badList=[]
    if bookDic.has_key(sMap.getName()) and bookDic[sMap.getName()].has_key("blacklist"):
        badSamples = bookDic[sMap.getName()]["blacklist"]
        badList= badSamples[:]
        for badSample in badSamples:
            des = sMap.getDescendants(badSample)
            for sample in des:
                if sample not in badSamples:
                    badList.append(sample)
    return badList
    
def cpProbeMaps(REALRUN,outDir,bookDic, sMap):
    if REALRUN in [-1,-2]:
        return
    
    dataPackageDir = outDir + sampleMapBaseName(sMap)
    for name in collectProbeMapBelongToSampleMap(bookDic, sMap.getName()):
        J = bookDic[name]
        if J['type'] != "probeMap":
            return False
        path = J['path']
        os.system("cp "+path+".json "+dataPackageDir+"/")
        if REALRUN:
            os.system("cp "+path+" "+dataPackageDir+"/")
    return

def sampleMapBaseName(sMap):
    baseName = sMap.getName()
    baseName = string.replace(baseName,".sampleMap","")
    baseName = string.replace(baseName,"_sampleMap","")
    baseName = string.replace(baseName,"sampleMap","")
    baseName = string.replace(baseName,"TCGA.","")
    return baseName
