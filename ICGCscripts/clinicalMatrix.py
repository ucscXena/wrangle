"""
Process a clinicalMatrix file.
"""
import sys, os, string, csv, json
from icgcLib import *
sys.path.insert(0,"../CGDataNew")
from ClinicalMatrixNew import *

def donor (fileIn, fileOut, xenaIds_specimen, xenaIds_donor, stat=None):
    print fileIn, fileOut

    fin = open(fileIn,'rU')
    fout = open(fileOut,'w')

    header= string.split(string.strip(fin.readline()),'\t')
    header[0]="xena_id"
    #header
    fout.write(header[0])
    for i in range(1,len(header)):
        fout.write("\t"+header[i])
    fout.write("\n")

    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        donor_id = data[0]
        xena_ids = xenaIds_donor[donor_id]
        for xena_id in xena_ids:
            fout.write(xena_id)
            for i in range(1,len(header)):
                fout.write("\t"+data[i])
            fout.write("\n")
    fin.close()
    fout.close()
    sortFileByFirstCol (fileOut)

    #survial
    if string.find(os.path.basename(fileIn),"donor.")!=-1:
        proj = findCohort(fileIn)
        survivalFile = os.path.dirname(fileOut)+"/survival."+proj+".tsv"
        survival(fileOut, survivalFile)
    return{}

def sortFileByFirstCol (fileIn):
    os.system("cp "+ fileIn +" .tmp")
    os.system("head -n 1 .tmp > "+fileIn)
    os.system("tail -n +2 .tmp | sort -t$'\t' >> " + fileIn)

def survival (fileIn, fileOut):
    print fileIn, fileOut

    fin = open(fileIn,'rU')

    col_donor_survival_time = -1
    col_donor_vital_status =-1 

    header= string.split(string.strip(fin.readline()),'\t')
    for i in range(1,len(header)):
        if header[i]=="donor_survival_time":
            col_donor_survival_time = i
        if header[i]=="donor_vital_status":
            col_donor_vital_status = i
    if col_donor_survival_time ==-1 or  col_donor_vital_status ==-1:
        os.system("rm -f "+ fileOut)
        return 

    fout = open(fileOut,'w')
    fout.write("xena_id\tdonor_survival_time\t_TIME_TO_EVENT\tdonor_vital_status\t_EVENT\n")
    
    has_data_event=False
    has_data_time=False

    for line in fin.readlines():
        data = string.split(line,'\t')
        id = data[0]
        time = data[col_donor_survival_time]
        event = data[col_donor_vital_status]
        if event =="alive":
            fout.write(string.join([id, time, time, event, "0"],"\t")+"\n")
        elif event =="deceased":
            fout.write(string.join([id, time, time, event, "1"],"\t")+"\n")
        else:
            fout.write(string.join([id, time, time, event, ""],"\t")+"\n")
        if time !="":
            has_data_time = True
        if event !="":
            has_data_event = True
    fin.close()
    fout.close()

    if has_data_time == False or has_data_event ==False:
        os.system("rm "+fileOut)
        print fileOut,"removed"
        return

    #json
    fout = open(fileOut+".json",'w')
    proj = findCohort(fileIn)
    xenaCohort = xenaCohortName(proj)
    obj= {
        ':clinicalFeature':"clinicalFeature",
        "version": date.today().isoformat(),
        "dataSubType":"phenotype",
        "type":"clinicalMatrix",
        'cohort': repInfo['name'] + ' ' + xenaCohort,
        'label':"Phenotype overall survival",
        'url': repInfo['download'] + repInfo['release'] + '/Projects/' + proj + '/donor.'+ proj+'.tsv.gz',
        'wrangling_procedure':"Donor data downloaded from dcc.icgc.org, extracted <b>donor_survival_time </b>and <b>donor_vital_status </b>information, loaded to UCSC xena database. donor_survival_time is _TIME_TO_EVENT.  donor_vital_status is _EVENT (deceased =1, alive =0).",
    }
    fout.write(json.dumps(obj, indent=4, sort_keys=True) + '\n')
    fout.close()
    sortFileByFirstCol (fileOut)
    return

def specimen (fileIn, fileOut, xenaIds_specimen, xenaIds_donor,  stat=None):
    print fileIn, fileOut

    fin = open(fileIn,'r')
    fout = open(fileOut,'w')

    ignoreList=[]
    dic={}
    header= string.split(string.strip(fin.readline()),'\t')
    header[0]="xena_id"
    for i in range(1,len(header)):
        if string.find(header[i],'specimen_id')!=-1:
            ignoreList.append(i)

    #header
    fout.write(header[0])
    for i in range(1,len(header)):
        if i not in ignoreList:
            fout.write("\t"+header[i])
    fout.write("\n")

    for line in fin.readlines():
        data = string.split(line,'\t')
        xena_id = xenaIds_specimen [data[0]]
        if xena_id not in dic:
            dic[xena_id]={}
        for i in range(1,len(header)):
            if i not in ignoreList:
                value = string.strip(data[i])
                if value in ["NA","na", "NULL","null","N/A"]:
                    value =""
                if header[i] not in dic[xena_id]:
                    dic[xena_id][header[i]]=[value]
                else:
                    if value not in dic[xena_id][header[i]]:
                        dic[xena_id][header[i]].append(value)

    for xena_id in dic.keys():
        fout.write(xena_id)
        for i in range(1, len(header)):
            if i not in ignoreList:
                column = header[i]
                if column in dic[xena_id]:
                    if len(dic[xena_id][column]) ==1:
                        fout.write("\t"+dic[xena_id][column][0])
                    else: # conflict remove
                        fout.write("\t")
                else:
                    fout.write("\t")
        fout.write("\n")
    fin.close()
    fout.close()
    sortFileByFirstCol (fileOut)

    return{}
