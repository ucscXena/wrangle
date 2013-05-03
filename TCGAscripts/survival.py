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

    #identify empty features
    badFeatures= finalClinMatrix.badCols()
    print "remove features", badFeatures
        
    if tag!=cancer:
        survivalMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_survival_"+tag)
    else:
        survivalMatrix= ClinicalMatrixNew(None,"clinical_"+cancer+"_survival")

    rO= overallSurvival (dir, finalClinMatrix, survivalMatrix, tag, cancer)
    rR = RFS (dir, finalClinMatrix, survivalMatrix, tag, cancer)
    if rO or rR:
        output (dir, finalClinMatrix,survivalMatrix, tag, cancer)
    
def overallSurvival (dir, finalClinMatrix, survivalMatrix, tag, cancer):
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
            if foundL:# and foundL!="0":
                survivalMatrix.setDATA(id,"_EVENT","0")
                survivalMatrix.setDATA(id,"_TIME_TO_EVENT",foundL)
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL_IND","0")
                #survivalMatrix.setDATA(id,"_OVERALL_SURVIVAL",foundL)
    return 1

def RFS  (dir, finalClinMatrix, survivalMatrix, tag, cancer):
    cols = finalClinMatrix.getCOLs()
    found=0
    for col in cols:
        if col=="person_neoplasm_cancer_status":
            found =1
            break
    if not found:
        print "no RFS info, can not compute _RFS, _RFS_ind"
        return 0
    
    foundDeath=0
    for col in cols:
        if col=="days_to_death":
            foundDeath =1
            break
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

    survivalMatrix.addOneColWithSameValue("_RFS","")
    survivalMatrix.addOneColWithSameValue("_RFS_IND","")

    minGood=0
    for id in finalClinMatrix.getROWs():
        #_EVENT         #_TIME_TO_EVENT
        foundD=0
        d = finalClinMatrix.getDATA(id, "days_to_new_tumor_event_after_initial_treatment")
        try:
            int(d)
            foundD =d 
            minGood= minGood+1
        except:
            # bad data no days_to_death for DECEASED
            pass
        d = finalClinMatrix.getDATA(id,"days_to_tumor_recurrence")
        try:
            int(d)
            if d > foundD:
                foundD=d
                minGood=minGood+1
        except:
            # bad data no days_to_death for DECEASED
            pass

        if foundD:
            survivalMatrix.setDATA(id,"_RFS_IND","1")
            survivalMatrix.setDATA(id,"_RFS",foundD)

        # if person_neoplasm_cancer_status = TUMOR FREE and no new tumor is detected above
        elif string.upper(finalClinMatrix.getDATA(id,"person_neoplasm_cancer_status")) =="TUMOR FREE":
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
            if foundDeath:
                d = finalClinMatrix.getDATA(id,"days_to_death") 
                try:
                    int(d)
                    if d>foundL:
                        foundL=d
                except:
                    pass
            if foundL:
                survivalMatrix.setDATA(id,"_RFS_IND","0")
                survivalMatrix.setDATA(id,"_RFS",foundL)
                
    if minGood<5:
        survivalMatrix.removeCols(["_RFS","_RFS_IND"])
    return 1
                
        
def output (dir, finalClinMatrix, survivalMatrix, tag, cancer):
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
    if tag!=cancer:
        cFfile =dir+"clinical_survival_clinicalFeature_"+tag
    else:
        cFfile =dir+"clinical_survival_clinicalFeature"
    fout = open(cFfile,"w")
    """
    fout.write("_OVERALL_SURVIVAL_IND\tshortTitle\toverall survivial indicator\n")
    fout.write("_OVERALL_SURVIVAL_IND\tlongTitle\t_OS_IND overall survivial indicator 0=censor (no_event) 1=event; derived from vital_status\n")
    fout.write("_OVERALL_SURVIVAL_IND\tvalueType\tcategory\n")

    fout.write("_OVERALL_SURVIVAL\tshortTitle\tOVERALL SURVIVAL\n")
    fout.write("_OVERALL_SURVIVAL\tlongTitle\t_OS overall survival; =days_to_death (if deceased); =max(days_to_last_known_alive,days_to_last_followup) (if living)\n")
    fout.write("_OVERALL_SURVIVAL\tvalueType\tfloat\n")
    """
    
    fout.write("_EVENT\tshortTitle\toverall survivial indicator\n")
    fout.write("_EVENT\tlongTitle\t_EVENT overall survivial indicator 0=censor (no_event) 1=event; derived from vital_status\n")
    fout.write("_EVENT\tvalueType\tcategory\n")

    fout.write("_TIME_TO_EVENT\tshortTitle\tOVERALL SURVIVAL\n")
    fout.write("_TIME_TO_EVENT\tlongTitle\t_TIME_TO_EVENT overall survival; =days_to_death (if deceased); =max(days_to_last_known_alive, days_to_last_followup) (if living)\n")
    fout.write("_TIME_TO_EVENT\tvalueType\tfloat\n")


    fout.write("_RFS_IND\tshortTitle\trecurrent free survival indicator\n")
    fout.write("_RFS_IND\tlongTitle\t_RFS_IND recurrent free survival indicator 1=new tumor; 0=otherwise and TUMOR FREE\n")
    fout.write("_RFS_IND\tvalueType\tcategory\n")

    fout.write("_RFS\tshortTitle\tRECURRENT FREE SURVIVAL\n")
    fout.write("_RFS\tlongTitle\t_RFS recurrent free survival; =max(days_to_new_tumor_event_after_initial_treatment, days_to_tumor_recurrence) (if event); =overall survival (if no event) \n")
    fout.write("_RFS\tvalueType\tfloat\n")

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

