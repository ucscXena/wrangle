"""
Process a clinicalMatrix file.
"""
import sys, os, string, csv, json
from icgcLib import *

cI = {
    'donor': 0,
    'specType': 1,
    'specId': 2,
    'sTime': 3,
    'sEvent': 4,
}

def initializeMetadata(fileIn):
    """
    Find the columns whose values are all empty
    and find the columns whose values are all the same.
    """
    with open(fileIn,'rU') as fIn:
        fIn = csv.DictReader(fIn, delimiter='\t')
        meta = fieldMetaInit(fIn.fieldnames, fIn.next())
        for row in fIn:
            fieldMetaRow(row, meta)
    meta['empty'] = filter(lambda f: meta['vals'][f] == '', meta['keys'])
    return meta

def buildHeader(meta):
    """
    Build the sort field order, without empty columns, putting these in front for
    sorting: donor, spec type, spec Id
    """
     # remove empty columns
    fields = filter(lambda f: not(f in meta['empty']), meta['all'])
    
    goodFields=[]
    # remove this because we combine all specIDs for each donorID + specType
    for field in fields:
        if string.find(field,'specimen_id')==-1:
            goodFields.append(field)

    fields= goodFields[:]
    # set up the header with constructed fields before the existing fields
    header = ['_SAMPLE_ID', '_PATIENT', 'icgc_specimen_id']
    if 'donor_survival_time' in fields and 'donor_vital_status' in fields:
        header += ['_TIME_TO_EVENT', '_EVENT']
        survivalHead = True
    else:
        survivalHead = False
    hFields = fields[:]
    header += hFields

    return {'header': header, 'fields': fields, 'survivalHead': survivalHead}

def nullDifferingVals(rows, prev): # TODO test
    """
    Within this donor and specimen type, null any values that differ between
    the entries.
    """
    indices = ['icgc_donor_id', 'specimenType', 'icgc_specimen_id']
    for col in prev.keys(): # prev is used to capture the information for this donor and specimen type
        if col in indices:
            continue

        val = rows[0][col] # initialize with first row of this column
        for row in rows[1:]:
            if row[col] != val:
                prev[col] = ''
                break

def buildXenaRow(prev, specIds, rows, xenaIds, fields, survivalHead):
    nullDifferingVals(rows, prev)
    xenaRow = [ xenaIds[prev['icgc_specimen_id']], prev['icgc_donor_id'], string.join(specIds, ', ') ]

    # these are combined into specIds, so drop their old fields
    del prev['icgc_specimen_id']
    newFields = fields[:]
    newFields.remove('icgc_specimen_id')

    # load our new survival values if both of the original survival fields have values
    if survivalHead:
        # if we ever want to clear our survival values if one is null...
        #if prev['donor_survival_time'] == '' or prev['donor_vital_status'] == '':
        #    xenaRow += ['', '']
        #else:
        #    xenaRow += [prev['donor_survival_time'], prev['donor_vital_status']]
        
        xenaRow += [prev['donor_survival_time'], prev['donor_vital_status']]

    return xenaRow + map(lambda x: prev[x], newFields)

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

    #survial
    if string.find(os.path.basename(fileIn),"donor.")!=-1:
        proj = findCohort(fileIn)
        survivalFile = os.path.dirname(fileOut)+"/survival."+proj+".tsv"
        survival(fileOut, survivalFile)

    return{}

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
        return 

    fout = open(fileOut,'w')
    fout.write("xena_id\tdonor_survival_time\t_TIME_TO_EVENT\tdonor_vital_status\t_EVENT\n")

    for line in fin.readlines():
        data = string.split(line,'\t')
        id = data[0]
        time = data[col_donor_survival_time]
        event = data[col_donor_vital_status]
        if event =="alive":
            fout.write(string.join([id, time, time, event, "0"],"\t")+"\n")
        elif event =="deceased":
            fout.write(string.join([id, time, time, event, "1"],"\t")+"\n")
    fin.close()
    fout.close()

    #json
    fout = open(fileOut+".json",'w')
    proj = findCohort(fileIn)
    xenaCohort = xenaCohortName(proj)
    obj= {
        "version": date.today().isoformat(),
        "dataSubType":"phenotype",
        "type":"clinicalMatrix",
        'cohort': repInfo['name'] + ' ' + xenaCohort,
        'label':"phenotype overall survival",
        ':clinicalFeature':"clinicalFeature",
    }
    fout.write(json.dumps(obj, indent=4, sort_keys=True) + '\n')
    fout.close()
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

    return{}
