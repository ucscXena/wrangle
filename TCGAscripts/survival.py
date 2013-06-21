import string, os, sys
import json,datetime
import csv

sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from CGDataLib import *
import  TCGAUtil

def survival (dir,cancer):
    os.system("rm "+dir+"clinical_survival")
    os.system("rm "+dir+"clinical_survival.json")
        
    ignore=0
    bookDic=cgWalk(dir,ignore)
    sampleMaps = collectSampleMaps(bookDic)
    missingMaps= collectMissingSampleMaps(bookDic)
    allMaps = sampleMaps.keys()
    allMaps.extend(missingMaps.keys())

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
        obj= bookDic[name]
        if obj['type']=="clinicalMatrix":
            #get matrix obj
            path = obj['path']
            name = obj['name']
            if string.find(name,"PANCAN")!=-1:
                continue
            cMatrix = ClinicalMatrixNew(path,name)
            
            if finalClinMatrix==None:
                finalClinMatrix= cMatrix
                
            #merge final and cMatrix
            if finalClinMatrix != cMatrix:
                r = finalClinMatrix.addNewCols(cMatrix,validation=True)
                if r!=True:
                    print "Fail to merge"
                    return False

    survivalMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_survival")

    rO= overallSurvival (dir, finalClinMatrix, survivalMatrix, cancer)
    rR = RFS (dir, finalClinMatrix, survivalMatrix, cancer)
    if rO or rR:
        output (dir, finalClinMatrix,survivalMatrix, cancer)
    
def overallSurvival (dir, finalClinMatrix, survivalMatrix, cancer):
    cols = finalClinMatrix.getCOLs()
    found =0
    for col in cols:
        if col=="vital_status":
            found =1
            break
    if not found:
        print "no vital_status, can not compute _TIME_TO_EVENT _EVENT"
        return 0
    found =0
    for col in cols:
        if col=="days_to_death":
            found =1
            break
    if not found:
        print "no days_to_death, can not compute _TIME_TO_EVENT _EVENT"
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

    survivalMatrix.addNewRows(finalClinMatrix.getROWs(),{})
    survivalMatrix.addOneColWithSameValue("_TIME_TO_EVENT","")
    survivalMatrix.addOneColWithSameValue("_EVENT","")
    #survivalMatrix.addOneColWithSameValue("_OVERALL_SURVIVAL","")
    #survivalMatrix.addOneColWithSameValue("_OVERALL_SURVIVAL_IND","")

    for id in finalClinMatrix.getROWs():
        #_EVENT         #_TIME_TO_EVENT
        v = finalClinMatrix.getDATA(id, "vital_status")
        if v=="DECEASED":
            d = finalClinMatrix.getDATA(id,"days_to_death")
            foundD =0
            try:
                int(d)
                foundD =1
                survivalMatrix.setDATA(id,"_EVENT","1")
                survivalMatrix.setDATA(id,"_TIME_TO_EVENT",d)
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL_IND","1")
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL",d)
                continue
            except:
                # bad data no days_to_death for DECEASED
                continue
        if v=="LIVING":
            if finalClinMatrix.getDATA(id,"FlagForSurvivalAnalysis") =="FLAG":
                print id,"FlagForSurvivalAnalysis"
                continue
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
                    if int(d)>int(foundL):
                        foundL=d
                except:
                    pass
            if foundL:# and foundL!="0":
                survivalMatrix.setDATA(id,"_EVENT","0")
                survivalMatrix.setDATA(id,"_TIME_TO_EVENT",foundL)
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL_IND","0")
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL",foundL)
    return 1

def RFS  (dir, finalClinMatrix, survivalMatrix, cancer):
    cols = finalClinMatrix.getCOLs()
    found=0
    for col in cols:
        if col=="new_tumor_event_after_initial_treatment":
            found =1
            break
    if not found:
        print "no new_tumor_event_after_initial_treatment info, can not compute _RFS, _RFS_ind"
        return 0
    
    survivalMatrix.addOneColWithSameValue("_RFS","")
    survivalMatrix.addOneColWithSameValue("_RFS_IND","")

    minGood=0
    for id in finalClinMatrix.getROWs():
        v = finalClinMatrix.getDATA(id, "new_tumor_event_after_initial_treatment")
        if v in ["YES","yes","Yes"]:
            foundD=0
            d = finalClinMatrix.getDATA(id,"days_to_new_tumor_event_after_initial_treatment")
            try:
                int(d)
                foundD = int(d)
            except:
                pass

            d = finalClinMatrix.getDATA(id,"days_to_tumor_recurrence")
            try:
                int(d)
                if int(d) > foundD:
                    foundD = int(d)
            except:
                pass

            if foundD >0:
                survivalMatrix.setDATA(id,"_RFS_IND","1")
                survivalMatrix.setDATA(id,"_RFS",str(foundD))
                minGood= minGood+1
                continue

        elif v in ["NO","no","No"]:
            d = survivalMatrix.getDATA(id,"_TIME_TO_EVENT") 
            try:
                int(d)
                survivalMatrix.setDATA(id,"_RFS_IND","0")
                survivalMatrix.setDATA(id,"_RFS",d)
            except:
                continue

    if minGood<5:
        survivalMatrix.removeCols(["_RFS","_RFS_IND"])
    return 1
                
        
def output (dir, finalClinMatrix, survivalMatrix, cancer):
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
    J["wrangling_procedure"]= "_EVENT from vital_status 0=no_event=LIVING 1=event=DECEASED; _TIME_TO_EVENT=days_to_death when _EVENT=1; _TIME_TO_EVENT=max(days_to_last_followup, days_to_last_known_alive) when _EVENT=0;  _RFS_IND 1=if there is days_to_new_tumor_event_after_initial_treatment or days_to_tumor_recurrenc information 0=otherwise and with person_neoplasm_status=TUMOR FREE; _RFS=max(days_to_new_tumor_event_after_initial_treatment, days_to_tumor_recurrence)"
    J["name"]=survivalMatrix.getName()
    J["type"]= "clinicalMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    
    #cFjson
    cFJ={}
    cFJ["name"]=J["name"]+"_clinFeat"
    cFJ["type"]="clinicalFeature"

    J[":clinicalFeature"] = cFJ["name"]

    #clinicalFeature
    cFfile =dir+"clinical_survival_clinicalFeature"
    fout = open(cFfile,"w")
    
    fout.write("_EVENT\tshortTitle\toverall survival indicator\n")
    fout.write("_EVENT\tlongTitle\t_EVENT overall survival indicator 0=censor (no_event) 1=event; derived from vital_status\n")
    fout.write("_EVENT\tvalueType\tcategory\n")

    fout.write("_TIME_TO_EVENT\tshortTitle\tOVERALL SURVIVAL\n")
    fout.write("_TIME_TO_EVENT\tlongTitle\t_TIME_TO_EVENT overall survival; =days_to_death (if deceased); =max(days_to_last_known_alive, days_to_last_followup) (if living)\n")
    fout.write("_TIME_TO_EVENT\tvalueType\tfloat\n")

    feature= "_EVENT"
    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(feature):
            priority= TCGAUtil.featurePriority[cancer][feature]
            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
            fout.write(feature+"\tvisibility\ton\n")

    feature= "_TIME_TO_EVENT"
    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(feature):
            priority= TCGAUtil.featurePriority[cancer][feature]
            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
            fout.write(feature+"\tvisibility\ton\n")

    fout.write("_RFS_IND\tshortTitle\trecurrence free survival indicator\n")
    fout.write("_RFS_IND\tlongTitle\t_RFS_IND recurrence free survival indicator 1=new tumor; 0=otherwise; derived from days_to_new_tumor_event_after_initial_treatment\n")
    fout.write("_RFS_IND\tvalueType\tcategory\n")

    fout.write("_RFS\tshortTitle\tRECURRENCE FREE SURVIVAL\n")
    fout.write("_RFS\tlongTitle\t_RFS recurrence free survival; =max(days_to_new_tumor_event_after_initial_treatment, days_to_tumor_recurrence) (if event); =overall survival (if no event) \n")
    fout.write("_RFS\tvalueType\tfloat\n")

    feature= "_RFS"
    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(feature):
            priority= TCGAUtil.featurePriority[cancer][feature]
            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
            fout.write(feature+"\tvisibility\ton\n")

    feature= "_RFS_IND"
    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(feature):
            priority= TCGAUtil.featurePriority[cancer][feature]
            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
            fout.write(feature+"\tvisibility\ton\n")
    fout.close()

    fout=open(dir+"clinical_survival.json",'w')
    fout.write( json.dumps( J, indent=-1 ) )
    fout.close()

    fout=open(dir+"clinical_survival_clinicalFeature.json",'w')
    fout.write( json.dumps( cFJ, indent=-1 ) )
    fout.close()

