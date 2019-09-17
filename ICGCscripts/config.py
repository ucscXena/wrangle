import urllib2,sys,json

#downloadDir = "/dev/shm/data/"
downloadDir = "/mnt/efsICGCnew/data/"
release = "release_28"

MAX_projects =100

icgcDataTypes = [ # only the dataset types of interest
    #'sample',
    #'specimen',
    #'donor',
    #'donor_exposure',
    #'donor_family',
    #'donor_therapy',

    #'protein_expression',  ## ALL US
    #'mirna_seq',
    #'meth_seq',
    #'meth_array',
    #'copy_number_somatic_mutation',
    #'simple_somatic_mutation.open',
    #'exp_seq',  ## NOT THE SAME UNIT
    #'structural_somatic_mutation',
    #'splice_variant'

    ##'exp_array',
]

def getProjects():
    projects=[]
    url = 'https://dcc.icgc.org/api/v1/projects?size='+ str(MAX_projects)
    response = urllib2.urlopen(url).read()
    J= json.loads(response)
    print len(J['hits'])
    for hit in J['hits']:
        projects.append(hit['id'])
    return projects

def getPrimarySite():
    projects = getProjects()
    dic={}
    for p in projects:
        url = 'https://dcc.icgc.org/api/v1/projects/'+p
        response = urllib2.urlopen(url).read()
        J= json.loads(response)
        try:
            dic [p]= J["primarySite"]
        except:
            dic[p]=""
    return dic

def getPrimaryDisease():
    projects = getProjects()
    dic={}
    for p in projects:
        url = 'https://dcc.icgc.org/api/v1/projects/'+p
        response = urllib2.urlopen(url).read()
        J= json.loads(response)
        try:
            disease = J["tumourType"]
        except:
            disease =""

        try:
            disease = disease +" : "+J["tumourSubtype"]
        except:
            pass
        dic[p]=disease
    return dic

