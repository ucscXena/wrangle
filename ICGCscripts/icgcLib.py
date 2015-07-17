"""
Utilities for wrangling, some specific to ICGC.
"""

import sys, os, gzip, string, csv, json, datetime, subprocess, smtplib
from itertools import *
from cfg import *
from repoMeta import *

projects = [ # ICGC projects
    'ALL-US',
    'AML-US',
    'BLCA-CN',
    'BOCA-FR',
    'BOCA-UK',
    'BRCA-EU',
    'BRCA-UK',
    'BTCA-SG',
    'CLLE-ES',
    'CCSK-US',
    'CMDI-UK',
    'COCA-CN',
    'EOPC-DE',
    'ESAD-UK',
    'ESCA-CN',
    'GACA-CN',
    'LAML-KR',
    'LIAD-FR',
    'LICA-FR',
    'LIHM-FR',
    'LINC-JP',
    'LIRI-JP',
    'LUSC-CN',
    'LUSC-KR',
    'MALY-DE',
    'MELA-AU',
    'NBL-US',
    'ORCA-IN',
    'OV-AU',
    'PACA-AU',
    'PACA-CA',
    'PACA-IT',
    'PAEN-AU',
    'PBCA-DE',
    'PRAD-CA',
    'PRAD-UK',
    'RECA-CN',
    'RECA-EU',
    'SKCA-BR',
    'THCA-SA',
    'WT-US'
]

icgcDataTypes = [ # only the dataset types of interest
    'specimen',
    'donor',
    'donor_exposure',
    'donor_family',
    'donor_therapy',
    'exp_array',
    'exp_seq',
    'mirna_seq',
    'meth_array',
    'simple_somatic_mutation.open',
    'mutGene'
]

dataSubTypes = { # referenced by repository data type name
    'donor':        {'name': 'phenotype',  'label':'phenotype donor', 'dsFx': 'donor',},
    'donor_family': {'name': 'phenotype',  'label':'phenotype family history', 'dsFx': 'donor',},
    'donor_therapy':{'name': 'phenotype',  'label':'phenotype therapy', 'dsFx': 'donor',},
    'donor_exposure':{'name': 'phenotype',  'label':'phenotype exposure','dsFx': 'donor',},
    'specimen':     {'name': 'phenotype',   'label':'phenotype specimen', 'dsFx': 'specimen',},
    'copy_number_somatic_mutation': {'name': 'copy number','label': 'copy number',             'dsFx': 'copy_number',},
    'exp_array':                    {'name': 'gene expression Array', 'label': 'gene expression Array',   'dsFx': 'gene_expression_array',},
    'exp_seq':                      {'name': 'gene expression RNAseq','label': 'gene expression RNAseq',  'dsFx': 'gene_expression_RNAseq',},
    'meth_array':                   {'name': 'DNA methylation','label': 'DNA methylation',         'dsFx': 'DNA_methylation',},
    'mirna_seq':                    {'name': 'miRNA expression RNAseq','label': 'miRNA expression RNAseq', 'dsFx': 'miRNA_expression_RNAseq',},
    'simple_somatic_mutation.open':      {'name': 'somatic mutation (SNPs and small INDELs)', 'label': 'somatic mutation (SNPs and small INDELs)',        'dsFx': 'mutationJson',},
    'mutGene':                      {'name': 'somatic non-silent mutation (gene-level)', 'label': 'somatic non-silent mutation (gene-level)', 'dsFx': 'mutGene',}
}

silent = False # Silent means silent execution requested, log only to file

repInfo  = { # ICGC repository information
    'release': 'release_19', 
    'version': date.today().isoformat(),
    'name': 'ICGC',
    'cohortDescrPrefix': 'This cohort is ',
    'dsDescrPrefix': 'This dataset is ',
    'download':'https://dcc.icgc.org/api/v1/download?fn=/',
}

repInfo['descrSuffix'] = ", "+ repInfo['release'] + '.'

def downloadUrlFromFile (repName, project, fileIn):
    if repName == "mutGene":
        outsideFile =string.replace(fileIn,'mutGene','simple_somatic_mutation.open')
    else:
        outsideFile = fileIn
    return repInfo['download'] + repInfo['release'] + '/Projects/' + project + '/' + outsideFile + '.gz'


def getCohortInfo():
    """
    Build cohort information from cohortMeta.tsv.
    """
    info = {}
    with open(rootDir + 'cohortMeta.tsv','rU') as fIn:
        fIn = csv.DictReader(fIn, delimiter='\t')

        fields = fIn.fieldnames
        for row in fIn:
            for field in fields:
                if field == 'repName':
                    repName = row[field]
                    info[repName] = {}
                if row[field] != None and row[field] != '':
                    info[repName][field] = row[field]
            info[repName]['metadataWritten'] = False
        #print json.dumps(info, indent=4, sort_keys=True)
        return info

def getDsInfo():
    """
    Build dataSubType information from dsMeta.tsv.
    """
    info = {}
    with open(rootDir + 'dsMeta.tsv','rU') as fIn:
        fIn = csv.DictReader(fIn, delimiter='\t')
        fields = fIn.fieldnames
        for row in fIn:
            for field in fields:
                if field == 'dataSubType':
                    dataSubType = row[field]
                    info[dataSubType] = {}
                elif row[field] != None and string.strip(row[field]) != '':
                    info[dataSubType][field] = row[field]
        return info

def initConfig():
    global cohortInfo
    cohortInfo = getCohortInfo()

def snpEff(eGene):
    try:
        return {'hGene': ensemblHugo[eGene], 'numFound': 1}
    except KeyError:
        return {'hGene': '', 'numFound': 0}

def loadSampleIdXref(proj):
    #xrefFile = dirs.ids + 'xref.' + proj
    xrefFile = dirs.ids + 'mapping.' + proj+".tsv"
    xref = {}
    with open(xrefFile,'rU') as fIn:
        fIn = csv.reader(fIn, delimiter='\t')
        for r in fIn:
            xref[r[0]] = r[1]
    return xref


def loadDonorIdXref(proj):
    #xrefFile = dirs.ids + 'xref.' + proj
    xrefFile = dirs.ids + 'mapping.' + proj+".tsv"
    xref = {}
    with open(xrefFile,'rU') as fIn:
        fIn = csv.reader(fIn, delimiter='\t')
        for r in fIn:
            if r[2] not in xref:
                xref[r[2]] = [r[1]]
            elif r[1] not in xref[r[2]]:
                xref[r[2]].append(r[1])
    return xref

def makeSilent(silentExec=True):
    global silent
    silent = silentExec

def isNumber(s):
    if s is None:
        return False
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def error(text):
    print 'ERROR: ' + text

def info(text):
    print 'INFO: ' + text

def infoTime(text):
    msg = str(datetime.datetime.now())[8:-7] + ' ' + text
    print 'INFO: ' + msg

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def writeClinicalFeature(cohort, dirOut):
    if cohortInfo[cohort]['metadataWritten']:
        return
    cohortInfo[cohort]['metadataWritten'] = True
    xenaCohort = xenaCohortName(cohort)
    meta = {
        'name': repInfo['name'] + ' ' + xenaCohort,
        'description': repInfo['cohortDescrPrefix'] + xenaCohort + repInfo['descrSuffix'],
        'version': repInfo['version'],
    }
    if 'primary_disease' in cohortInfo[cohort]:
        meta['primary_disease'] = cohortInfo[cohort]['primary_disease']
    fileOut = dirOut + '/cohort.' + cohort + '.json'
    with open(fileOut, 'w') as fOut:
        fOut.write(json.dumps(meta, indent=4, sort_keys=True) + '\n')

def getCohortName(cohort):
    return repInfo['name'] + '-' + cohortInfo[cohort]['repName']

def xenaCohortName(cohort):
    xenaName = cohortInfo[cohort]['xenaSuffix'] + ' (' + cohort + ')'
    return xenaName

def writeCohortMetadata(cohort, dirOut):
    if cohortInfo[cohort]['metadataWritten']:
        return
    cohortInfo[cohort]['metadataWritten'] = True
    xenaCohort = xenaCohortName(cohort)
    meta = {
        'name': repInfo['name'] + ' ' + xenaCohort,
        'description': repInfo['cohortDescrPrefix'] + xenaCohort + repInfo['descrSuffix'],
        'version': repInfo['version'],
    }
    if 'primary_disease' in cohortInfo[cohort]:
        meta['primary_disease'] = cohortInfo[cohort]['primary_disease']
    fileOut = dirOut + '/cohort.' + cohort + '.json'
    with open(fileOut, 'w') as fOut:
        fOut.write(json.dumps(meta, indent=4, sort_keys=True) + '\n')

def findCohort(fileIn):
    base = os.path.basename(fileIn)
    withoutTsv = base[:base.index('.tsv')]
    cohort = withoutTsv[withoutTsv.index('.')+1:]
    dotCount = cohort.count('.')
    if dotCount > 0:
        dotIdx = cohort.index('.')
        cohort = cohort[dotIdx+1:]
    return cohort

def findDataSubType(fileIn):
    base = os.path.basename(fileIn)
    reqDsName = base[:base.index('.')]

    if reqDsName == 'simple_somatic_mutation':
        reqDsName += '.open'
    try:
        dataSubTypeInfo = dataSubTypes[reqDsName]
    except KeyError:
        error('no dataSubType for file: ' + fileIn)
        return 'invalid'
    return dataSubTypeInfo

def findRepoName(fileIn):
    base = os.path.basename(fileIn)
    reqDsName = base[:base.index('.')]

    if reqDsName == 'simple_somatic_mutation':
        reqDsName += '.open'
    return reqDsName

def findRepDataType(dataSubType):
    for type in dataSubTypes:
        if dataSubTypes[type]['name'] == dataSubType:
            return type
    return None

def myArgParse(parser):
    parser.add_argument('--dirIn', metavar="Directory of original files from ICGC.", default=dirs.orig)
    parser.add_argument('--dirOut', metavar="Directory of xena-ready files.", default=dirs.xena)
    return parser.parse_args()

def buildDatasetCoreMetadata(fileIn, cohort):
    """
    Build the core dataset metadata given a dataSubType (dst) and cohort.
    """
    repName= findRepoName (fileIn)
    url = downloadUrlFromFile(repName, cohort, os.path.basename(fileIn))

    xenaCohort = xenaCohortName(cohort)

    repoInfo = repoMeta[repName] 
    label = repoInfo["label"]

    meta = {
        'cohort': repInfo['name'] + ' ' + xenaCohort,
        'description': repInfo['dsDescrPrefix'] + xenaCohort +  ' - ' + label + repInfo['descrSuffix'],
        'version': repInfo['version'],
        'url': url, 
    }

    LOG, value = alreadyLogged(fileIn)
    
    if repName == "simple_somatic_mutation.open":
        wrangle = "Data downloaded from dcc.icgc.org, converted to xena mutationVectorformats, loaded to UCSC xena database."

    elif repName == "mutGene":
        wrangle = "Data downloaded from dcc.icgc.org, converted to binary data to non-silent and non-silent mutations, binary results are loaded to UCSC xena database."

    elif repName in ["exp_array","exp_seq","mirna_seq"]:    
        if value =="":
            value = "downloaded value"
        else:
            value = value + " transformed value"
        if LOG: # already log transformed
            wrangle = "Data downloaded from dcc.icgc.org, converted to tab-delimited spread-sheet form, loaded to UCSC xena database."
            meta['unit']= value
        else:
            wrangle = "Data downloaded from dcc.icgc.org, converted tab-delimited spread-sheet/matrix form, values are log2(x+1) transformed, loaded to UCSC xena database."
            meta['unit']= "log2(x+1) of " + value 
    else:
        wrangle = "Data downloaded from dcc.icgc.org, converted to tab-delimited spread-sheet form, loaded to UCSC xena database."

    meta['wrangling_procedure'] = wrangle

    for key in repoInfo:
        meta[key]=repoInfo[key]

    return meta

def alreadyLogged (fileIn):
    LOG= False
    value = ""
    repName= findRepoName (fileIn)

    originalFile = dirs.orig +os.path.basename(fileIn)
    if not os.path.exists(originalFile):
        return [LOG, value]
    fin = open(originalFile,'rU') 
    fields = csv.DictReader(fin, delimiter='\t').fieldnames

    if 'normalization_algorithm' in fields:
        pos = fields.index('normalization_algorithm')
        for row in fin:
            value = string.split(row,'\t')[pos]
            if 'log' in value or 'RMA' in value:
                LOG= True
            fin.close()
            break

    return [LOG, value]

def fieldMetaInit(keys, row, used=None, useIdx=False):
    """
    Initialize the repository-specific metadata for a dataset.
    'all' are all of the fields names in the dataset
    'keys' are the field names to be stored in the metadata
    'vals' are used to determine if a column has all the same values, in which
        case that key and value is stored in the metadata
    'used' are the field names used in making xena-ready data
    'empty' are field names which lack a value in every row, & ignored for 
        xena-ready data & metadata
    """
    if useIdx:
        vals = {k: row[i] for i, k in enumerate(keys)}
    else:
        vals = {k: row[k] for k in keys}
    meta = {
        'keys': keys,
        'vals': vals,
        'all': keys[:],
    }
    if used:
        meta['used'] = used
    return meta


def fieldMetaRow(row, meta, useIdx=False):
    """
    Remove any fields from the repository-specific metadata for a dataset 
    when the column has more than one value.
    """
    def updateByIndex(key):
        i = meta['all'].index(key)
        return row[i] == meta['vals'][key]

    if useIdx:
        meta['keys'] = filter(updateByIndex, meta['keys'])
    else:
        meta['keys'] = filter(lambda x: row[x] == meta['vals'][x], meta['keys'])
    return meta

def buildMetadata(fields):
    # add the log2 message if log2 was applied
    if 'log2' in fields and fields['log2'] == True:
        if 'description' not in fields:
            fields['description']=""
        fields['description'] += ' The function log2(v+1) has been appied to all values.'

    return fields

def writeMetadata(fields, fileOut):
    """
    Write out the metadata for one dataset.
    """
    metadata = buildMetadata(fields)
    with open(fileOut + '.json', 'w') as fOut:
        fOut.write(json.dumps(metadata, indent=4, sort_keys=True) + '\n')

def initIcgcLib():
    global cohortInfo
    initConfig()

initIcgcLib()

