import string, os, sys
import json,datetime
import csv

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from CGDataLib import *

def survival (dir,cancer,tag):
    if tag!= cancer:
        os.system("rm "+dir+"clinical_survival_"+tag)
        os.system("rm "+dir+"clinical_survival_"+tag+".json")
    else:
        os.system("rm "+dir+"clinical_survival")
        os.system("rm "+dir+"clinical_survival.json")
        
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
                
    for name in datasets:  
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            if obj.has_key('upToDate') and obj['upToDate'] in ["yes", "Yes","YES"]:
                datasetsOrdered.insert(0,name)

    for name in datasetsOrdered:  
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

    cols = finalClinMatrix.getCOLs()
    found =0
    for col in cols:
        if col=="vital_status":
            found =1
            break
    if not found:
        print "no vital_status, can not compute _SURVIVAL _CENSOR"
        return 0
    found =0
    for col in cols:
        if col=="days_to_death":
            found =1
            break
    if not found:
        print "no days_to_death, can not compute _SURVIVAL _CENSOR"
        return 0

    foundAlive=0
    for col in cols:
        if col=="days_to_last_known_alive":
            foundAlive =1
            print "found days_to_last_known_alive"
            break
    foundFollowup=0
    for col in cols:
        if col=="days_to_last_followup":
            foundFollowup =1
            print "found days_to_last_followup"
            break

    if tag!=cancer:
        survivalMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_survival_"+tag)
    else:
        survivalMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_survival")
    survivalMatrix.addNewRows(finalClinMatrix.getROWs(),{})
    survivalMatrix.addOneColWithSameValue("_SURVIVAL","")
    survivalMatrix.addOneColWithSameValue("_CENSOR","")

    for id in finalClinMatrix.getROWs():
        #_CENSOR         #_SURVIVAL
        v = finalClinMatrix.getDATA(id, "vital_status")
        if v=="DECEASED":
            d = finalClinMatrix.getDATA(id,"days_to_death")
            foundD =0
            try:
                int(d)
                foundD =1
                survivalMatrix.setDATA(id,"_CENSOR","1")
                survivalMatrix.setDATA(id,"_SURVIVAL",d)
                continue
            except:
                # bad data no days_to_death for DECEASED
                continue
        if v=="LIVING":
            foundL =0
            if foundAlive:
                d = finalClinMatrix.getDATA(id,"days_to_last_known_alive") 
                try:
                    int(d)
                    foundL =d
                except:
                    pass
            if foundFollowup:
                d = finalClinMatrix.getDATA(id,"days_to_last_followup") 
                try:
                    int(d)
                    if d>foundL:
                        foundL=d
                except:
                    pass
            if foundL:
                survivalMatrix.setDATA(id,"_CENSOR","0")
                survivalMatrix.setDATA(id,"_SURVIVAL",foundL)
    if tag!=cancer:
        fout=open(dir+"clinical_survival_"+tag,'w')
    else:
        fout=open(dir+"clinical_survival",'w')
    survivalMatrix.store(fout,True)
    fout.close()

    #json
    J={}
    J["cgDataVersion"]=1
    J["redistribution"]= True
    J["dataProducer"]= "UCSC"
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
    J["wrangling_procedure"]= "_CENSOR=vital_status 0=no-event=LIVING 1=event=DECEASED; _SURVIVAL=days_to_death when _CENSOR=1; _SURVIVAL=max(days_to_last_followup, days_to_last_known_alive) when _CENSOR=0"
    J["name"]=survivalMatrix.getName()
    J["type"]= "clinicalMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #cFjson
    cFJ={}
    cFJ["name"]=J["name"]+"_clinFeat"
    cFJ["type"]="clinicalFeature"

    J[":clinicalFeature"] = cFJ["name"]

    #clinicalFeature
    if tag!=cancer:
        cFfile =dir+"clinical_survival_clinicalFeature_"+tag
    else:
        cFfile =dir+"clinical_survival_clinicalFeature"
    fout = open(cFfile,"w")
    fout.write("_CENSOR\tshortTitle\t_CENSOR\n")
    fout.write("_CENSOR\tlongTitle\t_CENSOR (from vital_status)\n")
    fout.write("_CENSOR\tvalueType\tcategory\n")

    fout.write("_SURVIVAL\tshortTitle\t_SURVIVAL\n")
    fout.write("_SURVIVAL\tlongTitle\t_SURVIVAL=days_to_death(event); _SURVIVAL=max(days_to_last_known_alive,days_to_last_followup) (no_event)\n")
    fout.write("_SURVIVAL\tvalueType\tfloat\n")
    fout.close()

    if tag!=cancer:
        fout=open(dir+"clinical_survival_"+tag+".json",'w')
    else:
        fout=open(dir+"clinical_survival.json",'w')
    fout.write( json.dumps( J, indent=-1 ) )
    fout.close()

    if tag!=cancer:
        fout=open(dir+"clinical_survival_clinicalFeature_"+tag+".json",'w')
    else:
        fout=open(dir+"clinical_survival_clinicalFeature.json",'w')
    fout.write( json.dumps( cFJ, indent=-1 ) )
    fout.close()
