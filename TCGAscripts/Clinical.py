import string, os, sys
import json,datetime
import csv
import datetime

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from survival import *
from ClinicalFeatureNew import *

#https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/bcr/biotab/clin/
#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/bcr/biotab/clin/

tmpDir="tmpTry/"

def ClinicalPublicBioTab(inDir, outDir, cancer,flog,REALRUN):
    PATHPATTERN= ""
    Clinical(inDir, outDir, cancer,flog,PATHPATTERN,REALRUN)
    
def Clinical(inDir, outDir, cancer,flog,PATHPATTERN, REALRUN):
    garbage=[tmpDir]
    if os.path.exists( tmpDir ):
        os.system("rm -rf "+tmpDir+"*")
    else:
        os.system("mkdir "+tmpDir)

    os.system("wget -r -l1 -H -nd -P"+ tmpDir +" -erobots=off  https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"+ string.lower(cancer) +"/bcr/biotab/clin/")

    dataDir=tmpDir
            
    #make sure there is data
    if REALRUN and (dataDir =="" or not os.path.exists(dataDir)):
        cleanGarbage(garbage)
        return

    process(inDir, outDir, dataDir, cancer,flog,PATHPATTERN, cancer,REALRUN)
    if REALRUN:
        survival(outDir+cancer+"/", cancer)

    if cancer in ["COAD","READ","LUAD","LUSC"]:
        if cancer in ["COAD","READ"]:
            deriveCancer="COADREAD"
        if cancer in ["LUAD","LUSC"]:
            deriveCancer="LUNG"
        process(inDir, outDir, dataDir, deriveCancer,flog,PATHPATTERN,  cancer,REALRUN)
        if REALRUN:
            survival(outDir+deriveCancer+"/", deriveCancer)

    cleanGarbage(garbage)
    return


def process (inDir, outDir, dataDir, cancer,flog,PATHPATTERN,  originCancer,REALRUN):
    #print status
    print cancer, __name__
    
    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    #data processing
    currentFollowUpV=0.0
    for file in os.listdir(dataDir):
        for pattern in ["clinical_sample","clinical_patient","clinical_follow_up","auxiliary","biospecimen_tumor_sample","biospecimen_sample"]:
            if string.find(file,pattern)!=-1:
                followUpV =0.0
                cgFileName = string.replace(file,".txt","")
                # the follow_up files has -vn.n version number
                if cgFileName!= re.sub(r'_v[1-9]+.[0-9]+','',cgFileName):
                    followUpV = float(string.split(cgFileName,"_")[3][1:])

                # the auxillary, biospecimen_tumor_sample file does not start with clin
                if cgFileName[0:9] != "clinical_":
                    cgFileName="clinical_"+cgFileName

                outfile = outDir+cancer+"/"+cgFileName
                cFfile =outfile+"_clinicalFeature"

                if not REALRUN:
                    if os.path.exists (cFfile):
                        tmpClinFeature= ClinicalFeatureNew(cFfile,"tmpName")
                        features= tmpClinFeature.getFeatures()
                        for feature in features:
                            if TCGAUtil.featurePriority.has_key(cancer):
                                if TCGAUtil.featurePriority[cancer].has_key(feature):
                                    priority= TCGAUtil.featurePriority[cancer][feature]
                                    tmpClinFeature.setFeaturePriority(feature, priority)
                                    tmpClinFeature.setFeatureVisibility(feature, "on")
                                    
                            stateOrder=None
                            if TCGAUtil.featureStateOrder.has_key(feature):
                                if TCGAUtil.featureStateOrder[feature].has_key(cancer):
                                    stateOrder = TCGAUtil.featureStateOrder[feature][cancer]
                                if TCGAUtil.featureStateOrder[feature].has_key("ALL"):
                                    stateOrder = TCGAUtil.featureStateOrder[feature]["ALL"]
                                print stateOrder
                            if stateOrder:
                                tmpClinFeature.setFeatureValueType(feature,"category")
                                tmpClinFeature.setFeatureStates(feature,stateOrder)
                                tmpClinFeature.setFeatureStateOrder(feature,stateOrder)
                                tmpClinFeature.setFeatureStateOrderRelax(feature,"true")

                        fout= open(cFfile,'w')    
                        tmpClinFeature.store(fout)
                        fout.close()
                        continue
                
                infile = dataDir+file
                #infile often row read has fewer fields than the fieldnames sequence
                # use csv.DictReader and Writer to fix this
                fin=open(infile,'r')
                reader = csv.DictReader(fin, delimiter="\t",restval="")
                fout=open(".tmp",'w')
                writer = csv.DictWriter(fout,delimiter="\t",fieldnames=reader.fieldnames)
                writer.writer.writerow(reader.fieldnames)
                writer.writerows(reader)
                fout.close()
                fin.close()
                os.system("cp .tmp "+infile)
                print pattern
                if pattern=="clinical_follow_up":
                    print file
                    if cancer ==originCancer:
                        cleanupFollowUpFile(infile, ".tmp")
                        os.system("cp .tmp "+infile)

                #clinicalMatrix
                if pattern =="biospecimen_tumor_sample":
                    FirstColAuto = True
                else:
                    FirstColAuto = False
                clinMatrix = ClinicalMatrixNew(infile, "foo", FirstColAuto)

                clinMatrix.removeCols(["ethnicity","race","jewish_origin","patient_id"])

                if pattern =="biospecimen_tumor_sample":
                    clinMatrix.removeCols(["vial_number"])
                    #replace bcr_sample_uuid with barcode
                    mapDic= TCGAUtil.uuid_Sample_all()
                    uuid_2_barcode (clinMatrix,"bcr_sample_uuid",mapDic,flog)
                    
                if pattern =="clinical_sample" or pattern=="biospecimen_sample":
                    if "sample_type" in clinMatrix.getCOLs():
                        add_col_PseudoSample (clinMatrix,"sample_type")
                    if "sample_type_id" in clinMatrix.getCOLs():
                        add_col_PseudoSample (clinMatrix,"sample_type_id")
                        
                #remove all cols with uuid
                features =clinMatrix.getCOLs()
                for f in features:
                    if string.find(f,"uuid")!=-1 or string.find(f,"UUID")!=-1 or string.find(f,"day_of")!=-1:
                        clinMatrix.removeCols([f])
                clinMatrix.replaceValue("null","")
                clinMatrix.replaceValue("NULL","")
                clinMatrix.replaceValue("Null","")
                clinMatrix.replaceValue("NA","")
                clinMatrix.replaceValue("[null]","")
                clinMatrix.replaceValue("[NULL]","")
                clinMatrix.replaceValue("[Null]","")
                clinMatrix.replaceValue("[NA]","")
                clinMatrix.replaceValue("[Not Available]|","")
                clinMatrix.replaceValue("[Not Available]","")
                clinMatrix.replaceValue("[Not Reported]","")
                clinMatrix.replaceValue("[Not Applicable]|","")
                clinMatrix.replaceValue("[Not Applicable]","")
                clinMatrix.replaceValue("[Not Requested]","")
                clinMatrix.replaceValue("[Completed]","")
                clinMatrix.replaceValue("[Pending]","")
                clinMatrix.replaceValue("Not Tested","")
                clinMatrix.replaceValue("[Not Evaluated]","")
                clinMatrix.replaceValue("[Unknown]","")
                clinMatrix.replaceValue("[]","")
                clinMatrix.replaceValue("LUNG","Lung") #stupid BCR
                clinMatrix.replaceValue("MSS|MSS","MSS") #stupid BCR
                clinMatrix.replaceValue("Alive","LIVING") #stupid BCR
                clinMatrix.replaceValue("ALIVE","LIVING") #stupid BCR
                clinMatrix.replaceValue("alive","LIVING") #stupid BCR
                clinMatrix.replaceValue("Dead","DECEASED") #stupid BCR
                clinMatrix.replaceValue("DEAD","DECEASED") #stupid BCR
                clinMatrix.replaceValue("dead","DECEASED") #stupid BCR

                #if cancer != originCancer:
                #    clinMatrix.addOneColWithSameValue("cohort",originCancer)
                    
                oHandle = open(outfile,"w")
                if pattern =="biospecimen_tumor_sample":
                    clinMatrix.storeSkip1stCol(oHandle, validation=True)
                else:
                    clinMatrix.store(oHandle, validation=True)
                oHandle.close()

                #clinicalFeature
                    
                fout = open(cFfile,"w")
                cFeatures = clinMatrix.getCOLs()
                for feature in cFeatures:
                    if not TCGAUtil.featureLongTitle.has_key(feature):
                        longTitle = feature
                        shortTitle = feature
                        message= "Feature Not in dictionary"+"\t"+feature+"\t"+feature
                        flog.write(message+"\n")
                    else:
                        longTitle = TCGAUtil.featureLongTitle[feature]
                        if TCGAUtil.featureShortTitle.has_key(feature):
                            shortTitle =TCGAUtil.featureShortTitle[feature]
                        else:
                            shortTitle =TCGAUtil.featureLongTitle[feature]

                    fout.write(feature+"\tshortTitle\t"+shortTitle+"\n")
                    fout.write(feature+"\tlongTitle\t"+longTitle+"\n")
                    if string.find(feature, "uuid")!=-1 or string.find(feature,"UUID")!=-1:
                        fout.write(feature+"\tvisibility\toff\n")
                    if TCGAUtil.valueType.has_key(feature):
                        fout.write(feature+"\tvalueType\t"+TCGAUtil.valueType[feature]+"\n")
                    stateOrder=None
                    if TCGAUtil.featureStateOrder.has_key(feature):
                        if TCGAUtil.featureStateOrder[feature].has_key(cancer):
                            fout.write(feature+"\tvalueType\tcategory\n")
                            stateOrder = TCGAUtil.featureStateOrder[feature][cancer]
                        if TCGAUtil.featureStateOrder[feature].has_key("ALL"):
                            fout.write(feature+"\tvalueType\tcategory\n")
                            stateOrder = TCGAUtil.featureStateOrder[feature]["ALL"]
                        if stateOrder:
                            for state in stateOrder:
                                fout.write(feature+"\tstate\t"+state+"\n")
                            fout.write(feature+"\tstateOrder\t\""+string.join(stateOrder,"\",\"")+"\"\n")
                            #if feature=="sample_type":
                            fout.write(feature+"\tstateOrderRelax\ttrue\n")
                            
                    if TCGAUtil.featurePriority.has_key(cancer):
                        if TCGAUtil.featurePriority[cancer].has_key(feature):
                            priority= TCGAUtil.featurePriority[cancer][feature]
                            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
                            fout.write(feature+"\tvisibility\ton\n")
                            
                    if feature in ["gender","age_at_initial_pathologic_diagnosis","days_to_last_followup","days_to_last_known_alive","sample_type","mononucleotide_and_dinucleotide_marker_panel_analysis_status","tumor_necrosis_percent","tumor_nuclei_percent"]:
                        fout.write(feature+"\tvisibility\ton\n")
                fout.close()

                #json
                J={}
                cFJ={}

                oHandle = open(outfile+".json","w")
                #stable
                if pattern =="clinical_sample":
                    if cancer != originCancer:
                        suffix = "clinSample"+PATHPATTERN+originCancer
                    else:
                        suffix = "clinSample"+PATHPATTERN
                if pattern =="clinical_patient":
                    if cancer != originCancer:
                        suffix = "clinPatient"+PATHPATTERN+originCancer
                    else:
                        suffix = "clinPatient"+PATHPATTERN
                if pattern =="clinical_follow_up":
                    if cancer != originCancer:
                        suffix = cgFileName+originCancer
                    else:
                        suffix = cgFileName
                if pattern =="auxiliary":
                    if cancer != originCancer:
                        suffix = "clinAuxiliary"+PATHPATTERN+originCancer
                    else:
                        suffix = "clinAuxiliary"+PATHPATTERN
                if pattern =="biospecimen_tumor_sample":
                    if cancer != originCancer:
                        suffix = "clinTumorSample"+PATHPATTERN+originCancer
                    else:
                        suffix = "clinTumorSample"+PATHPATTERN
                if  pattern=="biospecimen_sample":
                    if cancer != originCancer:
                        suffix = "bioSample"+PATHPATTERN+originCancer
                    else:
                        suffix = "bioSample"+PATHPATTERN
                J["cgDataVersion"]=1
                J["redistribution"]= True
                J["dataProducer"]= "TCGA biospecimen core resource"
                J["url"]=TCGAUtil.remoteBase \
                          +string.replace(inDir,TCGAUtil.localBase,"") \
                          + string.replace(dataDir,tmpDir,"")[:-1]
                J["version"]= datetime.date.today().isoformat()
                J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

                #change description
                J["wrangling_procedure"]= "Clinical data download from TCGA DCC, processed at UCSC into cgData repository"
                J["description"]= "This dataset is the TCGA "+ TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+") clinical data."
                
                #change cgData 
                J["name"]="TCGA_"+cancer+"_"+suffix

                cFJ["name"]=J["name"]+"_clinFeat"
                    
                cFJ["type"]="clinicalFeature"
                J["type"]= "clinicalMatrix" 
                J[":sampleMap"]="TCGA."+cancer+".sampleMap"
                J[":clinicalFeature"] = cFJ["name"]
                if pattern =="clinical_follow_up":
                    if cancer != originCancer:
                        J["upToDate"]=str(followUpV)+"_"+ originCancer  #"Yes"
                    else:
                        J["upToDate"]=str(followUpV)  #"Yes"
                if pattern =="biospecimen_tumor_sample":
                    J["outOfDate"]="Yes"
                oHandle.write( json.dumps( J, indent=-1 ) )
                oHandle.close()

                oHandle = open(cFfile+".json","w")
                oHandle.write( json.dumps( cFJ, indent=-1 ) )
                oHandle.close()
    return

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf "+dir+"*")
    return

def cleanupFollowUpFile(infile, outfile):
    fin = open(infile,'U')
    patients={}
    header =string.split(string.strip(fin.readline()),"\t")
    col_date_of_form_completion=-1
    col_days_to_last_followup =-1
    col_days_to_new_tumor_event_after_initial_treatment =-1
    flagPatientsForSurvival=[]
    for i in range (0,len(header)):
        if header[i]=="date_of_form_completion":
            col_date_of_form_completion=i
        if header[i]=="days_to_last_followup":
            col_days_to_last_followup =i
        if header[i]=="days_to_new_tumor_event_after_initial_treatment":
            col_days_to_new_tumor_event_after_initial_treatment=i
            
    for line in fin.readlines():
        data = string.split(line,'\t')
        sample= data[0]
        parts= string.split(sample,"-")
        patient = string.join(parts[0:3],"-")
        if len(parts)==4:
            findex=int(parts[3][1:])
        else:
            findex=0
        if not patients.has_key(patient):
            patients[patient]=[findex,data]
        elif findex > patients[patient][0]:
            #check sanity that the followup data makes sense, if not , set flag
            if col_date_of_form_completion != -1 and col_days_to_last_followup !=-1:
                old_date_of_form_completion = patients[patient][1][col_date_of_form_completion]
                new_date_of_form_completion = data[col_date_of_form_completion]
                delta_form = days_delta  (old_date_of_form_completion,new_date_of_form_completion)

                old_days_to_last_followup = patients[patient][1][col_days_to_last_followup]
                new_days_to_last_followup = data[col_days_to_last_followup]

                try:
                    delta_followup = int(new_days_to_last_followup) - int(old_days_to_last_followup)
                    if abs(delta_followup-delta_form) > 365:
                        #the patient need to be flagged for survival analysis, suspicious dates on followup
                        flagPatientsForSurvival.append(patient)
                except:
                    pass #contiue

            if col_days_to_new_tumor_event_after_initial_treatment !=-1:
                old_days_to_new_tumor_event_after_initial_treatment = patients[patient][1][col_days_to_new_tumor_event_after_initial_treatment]
                new_days_to_new_tumor_event_after_initial_treatment = data[col_days_to_new_tumor_event_after_initial_treatment]
                try:
                    old_days_to_new_tumor_event_after_initial_treatment = int(old_days_to_new_tumor_event_after_initial_treatment)
                except:
                    patients[patient]=[findex,data]
                    continue

                try:
                    new_days_to_new_tumor_event_after_initial_treatment = int(new_days_to_new_tumor_event_after_initial_treatment)
                    data[col_days_to_new_tumor_event_after_initial_treatment] = str(min(old_days_to_new_tumor_event_after_initial_treatment,new_days_to_new_tumor_event_after_initial_treatment))
                    
                except:
                    data[col_days_to_new_tumor_event_after_initial_treatment] = str(old_days_to_new_tumor_event_after_initial_treatment)
                print data[col_days_to_new_tumor_event_after_initial_treatment]
                patients[patient]=[findex,data]

            
    fin.close()
    
    fin = open(infile,'U')
    fout= open(outfile,'w')
    fout.write(string.strip(fin.readline())+"\tFlagForSurvivalAnalysis\n" )
    for line in fin.readlines():
        data = string.split(string.strip(line),'\t')
        sample= data[0]
        parts= string.split(sample,"-")
        patient = string.join(parts[0:3],"-")
        if len(parts)==4:
            findex=int(parts[3][1:])
        else:
            findex=0
        if findex == patients[patient][0]:
            if patient not in flagPatientsForSurvival:
                fout.write(patient+"\t"+string.join(data[1:],"\t")+"\t\n")
            else:
                fout.write(patient+"\t"+string.join(data[1:],"\t")+"\tFLAG\n")
    fout.close()

def days_delta  (old_date_of_form_completion,new_date_of_form_completion):
    y,m,d =string.split(old_date_of_form_completion,"-")
    y = int(y)
    m = int (m)
    d=int (d)
    if d==0:
        d=1

    dold = datetime.date(y,m,d)
    
    y,m,d =string.split(new_date_of_form_completion,"-")
    y = int(y)
    m = int (m)
    d=int (d)
    if d==0:
        d=1
    dnew = datetime.date(y,m,d)
    
    days_delta = (dnew-dold).days
    return days_delta

def uuid_2_barcode (clinMatrix,uuidcol,mapDic,flog): #convert uuid to barcode, if uuid not found, remove the sample
    rows= clinMatrix.getROWs()
    removeSamples=[]
    for row in rows:
        uuid = clinMatrix.getDATA(row,uuidcol)
        if TCGAUtil.is_barcode(uuid)==True:
            continue
        try:
            barcode = mapDic[string.lower(uuid)]
            clinMatrix.replaceValueInCol(uuidcol, uuid, barcode)
        except KeyError:
            removeSamples.append(row)
            print uuid,"not found"
            flog.write(uuid+" not found\n")

    if len(removeSamples)>0:
        r = clinMatrix.removeRows(removeSamples, True)
        if not r:
            print "fail to validate"
    clinMatrix.replaceColName(uuidcol,"tcgaBarCode")

def add_col_PseudoSample (clinMatrix,col): # add sample type informatin to pseudo samples
    rows= clinMatrix.getROWs()
    for row in rows:
        st = clinMatrix.getDATA(row,col)
        if st !=None and st!="":
            #assuming sample ids are TCGA barcode
            integration_id = TCGAUtil.barcode_IntegrationId(row)
            if clinMatrix.hasRow(integration_id):
                clinMatrix.setDATA(integration_id,col,st)
            else:
                clinMatrix.addNewRows([integration_id],{col:st})

    r=clinMatrix.validate()
    if r==False:
        print "add pseudoSample clinical infor", col, "fail"
