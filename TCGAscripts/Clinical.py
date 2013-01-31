import string, os, sys
import json,datetime
import csv

PATHPATTERN ="clinical"

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *
from CGDataUtil import *
from survival import *

#/inside/depot/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/*/bcr/minbiotab/clin

def ClinicalPublicBioTab(inDir, outDir, cancer,flog,REALRUN):
    PATHPATTERN= ""
    Clinical(inDir, outDir, cancer,flog,PATHPATTERN)
    
def ClinicalPublic(inDir, outDir, cancer,flog,REALRUN):
    PATHPATTERN= "public"
    Clinical(inDir, outDir, cancer,flog,PATHPATTERN)
    
def ClinicalAll(inDir, outDir, cancer,flog,REALRUN):
    PATHPATTERN = "all"
    Clinical(inDir, outDir, cancer,flog,PATHPATTERN)
    
def Clinical(inDir, outDir, cancer,flog,PATHPATTERN):
    garbage=["tmptmp/"]
    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")

    #single file in dir mode, uncompress to files
    dataDir =""
    lastDate=""
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue
        
        #file date
        lastDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)
        
        #is tar.gz?, uncompress 
        if string.find(file,".tar.gz")!=-1:
            os.system("tar -xzf "+inDir+file +" -C tmptmp/") 
            dataDir ="tmptmp/"
            break

    #make sure there is data
    if dataDir =="" or not os.path.exists(dataDir):
        cleanGarbage(garbage)
        return

    process(inDir, outDir, dataDir, cancer,flog,PATHPATTERN, lastDate, cancer)
    survival(outDir+cancer+"/", cancer, cancer)

    if cancer in ["COAD","READ"]:
        deriveCancer="COADREAD"
        process(inDir, outDir, dataDir, deriveCancer,flog,PATHPATTERN, lastDate, cancer)
        survival(outDir+deriveCancer+"/", deriveCancer, cancer)
    cleanGarbage(garbage)
    return


def process (inDir, outDir, dataDir, cancer,flog,PATHPATTERN, lastDate, originCancer):
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
                newCgFileName = re.sub(r'_v[1-9]+.[0-9]+','',cgFileName)
                if newCgFileName != cgFileName:
                    followUpV = float(string.split(cgFileName,"_")[3][1:])
                    if followUpV > currentFollowUpV:
                        currentFollowUpV = followUpV
                    
                cgFileName= newCgFileName
                # the auxillary, biospecimen_tumor_sample file does not start with clin
                if cgFileName[0:9] != "clinical_":
                    cgFileName="clinical_"+cgFileName
                    
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
                    if currentFollowUpV == followUpV:
                        cleanupFollowUpFile(infile, ".tmp")
                        os.system("cp .tmp "+infile)
                    else:
                        os.system("rm "+infile)
                        print followUpV, currentFollowUpV
                        continue
                outfile = outDir+cancer+"/"+cgFileName

                #clinicalMatrix
                if pattern =="biospecimen_tumor_sample":
                    FirstColAuto = True
                else:
                    FirstColAuto = False
                clinMatrix = ClinicalMatrixNew(infile, "foo", FirstColAuto)

                clinMatrix.removeCols(["ethnicity","race","jewish_origin"])

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
                clinMatrix.replaceValue("[Not Available]","")
                clinMatrix.replaceValue("[Not Reported]","")
                clinMatrix.replaceValue("[Not Applicable]","")
                clinMatrix.replaceValue("[Not Requested]","")
                clinMatrix.replaceValue("[Completed]","")
                clinMatrix.replaceValue("[Pending]","")
                clinMatrix.replaceValue("[]","")
                clinMatrix.replaceValue("Not Tested","")
                oHandle = open(outfile,"w")
                if pattern =="biospecimen_tumor_sample":
                    clinMatrix.storeSkip1stCol(oHandle, validation=True)
                else:
                    clinMatrix.store(oHandle, validation=True)
                oHandle.close()

                #clinicalFeature
                cFfile =outfile+"_clinicalFeature"
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
                    if TCGAUtil.featureStateOrder.has_key(feature):
                        fout.write(feature+"\tvalueType\tcategory\n")
                        stateOrder = TCGAUtil.featureStateOrder[feature]                        
                        for state in stateOrder:
                            fout.write(feature+"\tstate\t"+state+"\n")
                        fout.write(feature+"\tstateOrder\t\""+string.join(stateOrder,"\",\"")+"\"\n")

                    if TCGAUtil.featurePriority.has_key(cancer):
                        if TCGAUtil.featurePriority[cancer].has_key(feature):
                            priority= TCGAUtil.featurePriority[cancer][feature]
                            fout.write(feature+"\tpriority\t"+str(priority)+"\n")
                            fout.write(feature+"\tvisibility\ton\n")
                            
                    if feature in ["age_at_initial_pathologic_diagnosis","days_to_last_followup","days_to_last_known_alive","sample_type","mononucleotide_and_dinucleotide_marker_panel_analysis_status","cancer type","tumor_necrosis_percent","tumor_nuclei_percent"]:
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
                        suffix = "clinFollowup"+PATHPATTERN+originCancer
                    else:
                        suffix = "clinFollowup"+PATHPATTERN
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
                          + string.replace(dataDir,"tmptmp/","")[:-1]
                J["version"]= datetime.date.today().isoformat()
                J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()

                #change description
                J["wrangling_procedure"]= "Clinical data download from TCGA DCC, processed at UCSC into cgData repository"
                J["description"]= "This dataset is the TCGA "+ TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+") clinical data."
                
                #change cgData 
                J["name"]="TCGA_"+cancer+"_"+suffix
                name = trackName_fix(J['name'])
                if name ==False:
                    message = "bad object name, need fix otherwise break loader, too long "+J["name"]
                    print message
                    flog.write(message+"\n")
                    return
                elif not trackName_good(name+"_clinFeat"):
                    message = "bad object name, need fix otherwise break loader, too long "+J["name"]+"-clinicalFeature"
                    print message
                    flog.write(message+"\n")
                    return
                else:
                    J["name"]=name

                cFJ["name"]=J["name"]+"_clinFeat"
                    
                cFJ["type"]="clinicalFeature"
                J["type"]= "clinicalMatrix" 
                J[":sampleMap"]="TCGA."+cancer+".sampleMap"
                J[":clinicalFeature"] = cFJ["name"]
                if pattern =="clinical_follow_up":
                    J["upToDate"]="Yes"
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
        os.system("rm -rf dir")
    return

def cleanupFollowUpFile(infile, outfile):
    fin = open(infile,'U')
    patients={}
    fin.readline()
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
            patients[patient]=findex
        elif findex > patients[patient]:
            patients[patient]=findex            
    fin.close()
    
    fin = open(infile,'U')
    fout= open(outfile,'w')
    fout.write(fin.readline())
    for line in fin.readlines():
        data = string.split(line,'\t')
        sample= data[0]
        parts= string.split(sample,"-")
        patient = string.join(parts[0:3],"-")
        if len(parts)==4:
            findex=int(parts[3][1:])
        else:
            findex=0
        if findex == patients[patient]:
            fout.write(patient+"\t"+string.join(data[1:],"\t"))
    fout.close()

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
