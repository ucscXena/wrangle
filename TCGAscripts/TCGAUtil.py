#TCGA Barcode to UUID Web Service User's Guide to allow Barcode-UUID lookup  https://wiki.nci.nih.gov/x/m4n9AQ
#UUID service https://wiki.nci.nih.gov/display/TCGA/TCGA+Biospecimen+Metadata+Web+Service+User's+Guide
#UUID browser https://tcga-data.nci.nih.gov/uuid/uuidBrowser.htm

#UUID Migration Plan Wiki page at: https://wiki.nci.nih.gov/x/xh9hAg
#TCGA Variant Call Format (VCF) 1.1 Specification https://wiki.nci.nih.gov/x/2gcYAw
#Sample and Data Relationship Format https://wiki.nci.nih.gov/x/9aFXAg
#Mutation Annotation Format (MAF) Specification https://wiki.nci.nih.gov/x/eJaPAQ

import string
import json
import os
import os.path, time,datetime

from TCGAfeature import *

remoteBase="https://tcga-data.nci.nih.gov/"
localBase="/inside/depot/"

uuid_barcode_txt= "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/metadata/metadata.20130228.txt"

# all sorted by short titles
cancerOfficial={
    "LAML":"acute myeloid leukemia",
    "BLCA":"bladder urothelial carcinoma",
    "LGG":"brain lower grade glioma",
    "BRCA":"breast invasive carcinoma",
    "CESC":"cervical/endocervical squamous cell carcinoma",
    "CNTL":"controls",
    "COAD":"colon adenocarcinoma",
    "COADREAD":"colon & rectum adenocarcinoma",
    "ESCA":"esophageal carcinoma",
    "GBM":"glioblastoma multiforme",
    "HNSC":"head & neck squamous cell carcinoma",
    "KICH":"kidney chromophobe",
    "KIRC":"kidney renal clear cell carcinoma",
    "KIRP":"kidney renal papillary cell carcinoma",
    "LIHC":"liver hepatocellular carcinoma",
    "LUAD":"lung adenocarcinoma",
    "LUSC":"lung squamous cell carcinoma",
    "DLBC":"lymphoid neoplasm diffuse large B-cell lymphoma",
    "LNNH":"lymphoid neoplasm non-Hodgkins lymphoma",
    "OV":"ovarian serous cystadenocarcinoma",
    "PAAD":"pancreatic adenocarcinoma",
    "PRAD":"prostate adenocarcinoma",
    "READ":"rectum adenocarcinoma",
    "SARC":"sarcoma",
    "SKCM":"skin cutaneous melanoma",
    "STAD":"stomach adenocarcinoma",
    "THCA":"thyroid carcinoma",
    "UCEC":"uterine corpus endometrioid carcinoma",
    "ACC":"adrenocortical carcinoma",
    "LCLL":"chronic lymphocytic leukemia",
    "LCML":"chronic myelogenous leukemia",
    "MESO":"mesothelioma",
    "MISC":"miscellaneous",
    "PCPG":"pheochromocytoma and paraganglioma",
    "UCS":"uterine carcinosarcoma"
    }

cancerHumanReadable={
    "LAML":"Acute Myeloid Leukemia",
    "BLCA":"Bladder Cancer",
    "LGG":"Lower Grade Glioma",
    "BRCA":"Breast Cancer",
    "CESC":"Cervical Cancer",
    "CNTL":"Controls",
    "COAD":"Colon Cancer",
    "COADREAD":"Colon and Rectal Cancer",
    "ESCA":"Esophageal Caner",
    "GBM":"Glioblastoma",
    "HNSC":"Head and Neck Cancer",
    "KICH":"Kidney Chromophobe",
    "KIRC":"Kidney Clear Cell Carcinoma",
    "KIRP":"Kidney Papillary Cell Carcinoma",
    "LIHC":"Liver Cancer",
    "LUAD":"Lung Adenocarcinoma",
    "LUSC":"Lung Squamous Cell Carcinoma",
    "DLBC":"Large B-cell Lymphoma",
    "LNNH":"Non-Hodgkins Lymphoma",
    "OV":"Ovarian Cancer",
    "PAAD":"Pancreatic Cancer",
    "PRAD":"Prostate Cancer",
    "READ":"Rectal Cancer",
    "SARC":"Sarcoma",
    "SKCM":"Melanoma",
    "STAD":"Stomach Cancer",
    "THCA":"Thyroid Cancer",
    "UCEC":"Endometrioid Cancer",
    "ACC":"Adrenocortical Cancer",
    "LCLL":"Chronic Lymphocytic Leukemia",
    "LCML":"Chronic Myeloid Leukemia",
    "MESO":"Mesothelioma",
    "MISC":"Miscellaneous",
    "PCPG":"Pheochromocytoma & Paraganglioma",
    "UCS":"Uterine Cancer"
    }

cancerGroupTitle={
    "LAML":"acute myeloid leukemia",
    "BLCA":"bladder urothelial carcinoma",
    "LGG":"brain lower grade glioma",
    "BRCA":"breast invasive carcinoma",
    "CESC":"cervical and endocervical SCC",
    "CNTL":"controls",
    "COAD":"colon adenocarcinoma",
    "COADREAD":"colon & rectum adenocarcinoma",
    "ESCA":"esophageal carcinoma",
    "GBM":"glioblastoma multiforme",
    "HNSC":"head & neck squamous cell carcinoma",
    "KICH":"kidney chromophobe",
    "KIRC":"kiney clear cell carcinoma",
    "KIRP":"kiney papillary cell carcinoma",
    "LIHC":"liver hepatocellular carcinoma",
    "LUAD":"lung adenocarcinoma",
    "LUSC":"lung squamous cell carcinoma",
    "DLBC":"diffuse large B-cell lymphoma",
    "LNNH":"non-Hodgkins lymphoma",
    "OV":"ovarian serous cystadenocarcinoma",
    "PAAD":"pancreatic adenocarcinoma",
    "PRAD":"prostate adenocarcinoma",
    "READ":"rectum adenocarcinoma",
    "SARC":"sarcoma",
    "SKCM":"skin cutaneous melanoma",
    "STAD":"stomach adenocarcinoma",
    "THCA":"thyroid carcinoma",
    "UCEC":"uterine corpus endometrioid carcinoma",
    "ACC":"adrenocortical cancer",
    "LCLL":"chronic lymphocytic leukemia",
    "LCML":"chronic myeloid leukemia",
    "MESO":"mesothelioma",
    "MISC":"miscellaneous",
    "PCPG":"pheochromocytoma & paraganglioma",
    "UCS":"uterine Cancer"
    }

TCGASampleType={
    "01":"Primary solid Tumor",
    "02":"Recurrent Solid Tumor",
    "03":"Primary Blood Derived Cancer - Peripheral Blood",
    "04":"Recurrent Blood Derived Cancer - Bone Marrow",
    "05":"Additional - New Primary",
    "06":"Metastatic",
    "07":"Additional Metastatic",
    "08":"Human Tumor Original Cells",
    "09":"Primary Blood Derived Cancer - Bone Marrow",
    "10":"Blood Derived Normal",
    "11":"Solid Tissue Normal",
    "12":"Buccal Cell Normal",
    "13":"EBV Immortalized Normal",
    "14":"Bone Marrow Normal",
    "20":"Control Analyte",
    "40":"Recurrent Blood Derived Cancer - Peripheral Blood",
    "50":"Cell Lines",
    "60":"Primary Xenograft Tissue",
    "61":"Cell Line Derived Xenograft Tissue"
}

valueType={
    "patient_id":"category",
    "":""
    }

featureStateOrder={
    "er_level_cell_percentage_category":{"BRCA":
                                         ["<10%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-99%"]},
    "her2_erbb_pos_finding_cell_percent_category":{"BRCA":
                                                   ["<10%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-99%"]},
    "distant_metastasis_pathologic_spread":{"ALL":
                                            ['MX', 'M0', 'M1', 'M1a', 'M1b']},
    "histological_type":{"LGG":
                         ["Oligodendroglioma","Oligoastrocytoma","Astrocytoma"]},
    "neoplasm_histologic_grade":{"ALL":
                                 ["GX","GB","G1","G2","G3","G4"]},
    "sample_type":{"ALL":
                   ["Primary Tumor","Primary solid Tumor",
                   "Primary Blood Derived Cancer - Peripheral Blood",
                   "Primary Blood Derived Cancer - Bone Marrow",
                   "Recurrent Tumor",
                   "Recurrent Solid Tumor",
                   "Recurrent Blood Derived Cancer - Peripheral Blood",
                   "Recurrent Blood Derived Cancer - Bone Marrow",
                   "Additional - New Primary",
                   "Metastatic",
                   "Additional Metastatic", 
                   "Human Tumor Original Cells", 
                   "Blood Derived Normal",
                   "Solid Tissue Normal",
                   "Buccal Cell Normal",
                   "EBV Immortalized Normal",
                   "Bone Marrow Normal",
                   "Cell Lines",
                   "Primary Xenograft Tissue",
                   "Cell Line Derived Xenograft Tissue",
                   "Control Analyte"]},
    "tobacco_smoking_history_indicator":{"LUSC":
                                         ["Lifelong Non-smoker",
                                         "Current reformed smoker for > 15 years",
                                         "Current reformed smoker for < or = 15 years",
                                         "Current smoker"],
                                         "LUAD":
                                         ["Lifelong Non-smoker",
                                         "Current reformed smoker for > 15 years",
                                         "Current reformed smoker for < or = 15 years",
                                         "Current smoker"]
                                         },
    "vital_status":{"ALL":
                    ["LIVING","DECEASED"]},
    "":""
    }

featurePriority={
    "BLCA": {"sample_type":"1",
             "":""    
             },
    "BRCA": {"sample_type":"1",
             "PAM50Call":"2",
             "er_level_cell_percentage_category":"3",
             "her2_immunohistochemistry_level_result":"4",
             "":""
             },
    "COAD": {"braf_gene_analysis_result":"1",
             "kras_mutation_found":"2",
             "CIMP":"3",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"4",
             "hypermutation":"5",
             "":""
             },
    "COADREAD":{"tissue_source_site":"1.0",
                "hypermutation":"1.1",
                "braf_gene_analysis_result":"1.2",
                "kras_mutation_found":"1.3",
                "CIMP":"1.4",
                "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"1.5",
                "":""
                },
    "CESC": {"sample_type":"1",
             "neoplasm_histologic_grade":"2",
             "":""
             },
    "ECSA": {"sample_type":"1",
             "":""
             },
    "DLBC": {"sample_type":"1",
             "":""
             },
    "GBM":  {"sample_type":"1",
             "GeneExp_Subtype":"2",
             "G_CIMP_STATUS":"3",
             "CDE_survival_time":"4",
             "":""    
             },
    "HNSC": {"sample_type":"1",
             "tumor_stage":"2",
             "age_at_initial_pathologic_diagnosis":"3",
             "days_to_last_known_alive":"4",
             "":""    
             },
    "KIRC": {"sample_type":"1",
             "":""    
             },
    "KIRP": {"sample_type":"1",
             "":""    
             },
    "LIHC": {"sample_type":"1",
             "neoplasm_histologic_grade":"2",
             "ajcc_tumor_stage_code":"3",
             "":""    
             },
    "LAML": {"sample_type_id":"1",
             "acute_myeloid_leukemia_calgb_cytogenetics_risk_category":"2",
             "gender":"3",
             "age_at_initial_pathologic_diagnosis":"4",
             "":""    
             },
    "LGG":  {"sample_type":"1",
             "histological_type":"2",
             "histologic_classification":"3",
             "days_to_last_followup":"4",
             "":""    
             },
    "LUAD": {"sample_type":"1",
             "Expression_Subtype":"1.2",
             "EGFR":"1.3",
             "KRAS":"1.4",
             "tobacco_smoking_history_indicator":"2",
             "tumor_stage":"3",
             "distant_metastasis_pathologic_spread":"4",
             "":""    
             },
    "LUSC": {"sample_type":"1",
             "":""    
             },
    "KICH": {"sample_type":"1",
             "":""    
             },
    "OV": {"sample_type":"1",
           "neoplasm_histologic_grade":"2",
           "":"" 
           },
    "PAAD": {"sample_type":"1",
           "":""    
           },
    "PRAD": {"sample_type":"1",
             "":""    
             },
    "READ": {"braf_gene_analysis_result":"1.0",
             "kras_mutation_found":"1.1",
             "CIMP":"1.2",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"1.3",
             "hypermutation":"1.4",
             "":""
             },
    "SKCM": {"sample_type":"1",
             "gender":"2",
             "":""    
             },
    "STAD": {"sample_type":"1",
             "":""    
             },
    "THCA": {"sample_type":"1",
             "":""    
             },
    "UCEC": {"sample_type":"1",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"2",
             "":""
             },
    "ACC":{"sample_type":"1",
             "":""    
             },
    "LCLL":{"sample_type":"1",
             "":""    
             },
    "LCML":{"sample_type":"1",
             "":""    
             },
    "MESO":{"sample_type":"1",
             "":""    
             },
    "MISC":{"sample_type":"1",
             "":""    
             },
    "PCPG":{"sample_type":"1",
             "":""    
             },
    "UCS":{"sample_type":"1",
             "":""    
             }
    }


clinDataDesc ="The accompanied clinical data are downloaded from TCGA data coordination center. They are the patient, sample, follow up, biospecimen tumor sample, and auxiliary clinical information. Names of the clinical fields were \"curated\" by UCSC to be more readable."

def is_barcode (id):
    if len(id)<12:
        return False
    if id[:4]!="TCGA":
        return False
    if id[4]!="-" and id[7]!="-":
        return False
    if string.split(id,"-") < 3:
        return False
    return True
    
def barcode_SampleType (barcode):
    #return the sample type code such as 01 (is primary tissue) and  20 (is cell line control)
    parts= string.split(barcode,"-")
    if len(parts)>=4:
        return parts[3][:2]
    else:
        return False

def barcode_IntegrationId (barcode):
    #return TCGA-01-1234
    parts= string.split(barcode,"-")
    if len(parts)<4:
        return None
    if len(parts[3])>2:
        parts[3]=parts[3][0:2]
    return string.join(parts[0:4],"-")
    
## UUID related code
def uuid_normal_cellline():
    """
    txtlocal =  "uuid_normal_cellline.txt"
    dic={}
    fin = open(txtlocal)
    for line in fin.readlines():
        dic[string.strip(line)]=0
    fin.close()
    return dic
    """

    txtlocal =  "uuid_normal_cellline.txt"
    localfile = "uuid_normal_cellline.json"
    
    current=0
    if os.path.exists(txtlocal):
        fileTime= time.strftime("%Y-%m-%d",time.localtime(os.path.getmtime(localfile)))
        today= datetime.date.today()
        if fileTime ==str(today):
            current=1
    if not current:
        site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/"
        q="json?sampleType=Normal,cellLine"  #Possible values: 'Tumor' (codes 01-09), 'Normal' (codes 10-19), 'cellLine' (code 20) 
        os.system("cp "+localfile +" "+localfile+"_BK")
        os.system("wget -nv -O "+localfile +" "+site+q)
        
        iHandle = open(localfile,"r")
        J = json.loads(iHandle.read())
        iHandle.close()

        dic={}
        for item in J["tcgaElement"]:
            address = item["@href"]
            uuid = string.split(address,"/")[-1]
            dic[uuid]=0
        #write to txtlocal
        os.system("cp "+txtlocal+" "+txtlocal+"_BK")
        fout =open(txtlocal,"w")
        for index in dic:
            fout.write(index+'\n')
        fout.close()
        return dic

    else: #parse texlocal
        dic={}
        fin = open(txtlocal)
        for line in fin.readlines():
            dic[string.strip(line)]=0
        fin.close()
        return dic

        
def uuid_cellline():
    """
    txtlocal = "uuid_cellline.txt"
    
    dic={}
    fin = open(txtlocal)
    for line in fin.readlines():
        dic[string.strip(line)]=0
    fin.close()
    return dic
    """

    txtlocal = "uuid_cellline.txt"
    localfile = "uuid_cellline.json"
    current=0
    if os.path.exists(txtlocal):
        fileTime= time.strftime("%Y-%m-%d",time.localtime(os.path.getmtime(localfile)))
        today= datetime.date.today()
        if fileTime ==str(today):
            current=1
    if not current:
        site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/json?sampleType=20"
        os.system("cp "+localfile+" "+localfile+"_BK")
        os.system("wget -nv -O "+localfile +" "+site)
        
        iHandle = open(localfile,"r")
        J = json.loads(iHandle.read())
        iHandle.close()

        dic={}
        for item in J["tcgaElement"]:
            address = item["@href"]
            uuid = string.split(address,"/")[-1]
            dic[uuid]=0

        #write to txtlocal
        os.system("cp "+txtlocal+" "+txtlocal+"_BK")
        fout =open(txtlocal,"w")
        for index in dic:
            fout.write(index+'\n')
        fout.close()
        return dic
    
    else: #parse texlocal
        dic={}
        fin = open(txtlocal)
        for line in fin.readlines():
            dic[string.strip(line)]=0
        fin.close()
        return dic
        
def uuid_Aliquot_all():  #the dictionary of all Aliquote key=uuid value=barcode
    localfile = "uuid_Aliquot.json"
    current=0
    if os.path.exists(localfile):
        fileTime= time.strftime("%Y-%m-%d",time.localtime(os.path.getmtime(localfile)))
        today= datetime.date.today()
        if fileTime ==str(today):
            current=1
    if not current:
        site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/"
        q="json?elementType=Aliquot"
        os.system("wget -nv -O "+localfile +" "+site+q)
        
    iHandle = open(localfile,"r")
    J = json.loads(iHandle.read())
    iHandle.close()

    dic={}
    for item in J["tcgaElement"]:
        address = item["@href"]
        uuid = string.split(address,"/")[-1]
        dic[uuid]=0

    uuid_barcode_dic={}
    os.system("wget -O uuid_barcode "+uuid_barcode_txt)
    fin = open("uuid_barcode","r")
    for line in fin.readlines():
        #uuid, barcode =string.split(line,',')[0:2]
        uuid, barcode =string.split(line,'\t')[0:2]
        if not uuid_barcode_dic.has_key(uuid):
            uuid_barcode_dic[uuid]=barcode
    fin.close()

    if os.path.exists("jing_uuid_barcode"):
        fin = open("jing_uuid_barcode","r")
        for line in fin.readlines():
            line =line[:-1]
            if line =="":
                continue
            uuid, barcode =string.split(line,'\t')[0:2]
            if not uuid_barcode_dic.has_key(uuid):
                uuid_barcode_dic[uuid]=barcode
            uuid_barcode_dic[uuid]=barcode
        fin.close()
        
    for item in dic:
        if item not in uuid_barcode_dic:
            c=1
            while c:
                try:
                    site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/"
                    q="json/uuid/"+item
                    os.system("wget -nv -O tmp "+ site+q)
                    iHandle = open("tmp","r")
                    J = json.loads(iHandle.read())
                    iHandle.close()
                    uuid_barcode_dic[item]=J["tcgaElement"]["barcodes"]["barcode"]
                    break
                except ValueError:
                    c=c+1
                    if c>3:
                        break

    os.system("cp jing_uuid_barcode jing_uuid_barcode_BK")
    fout= open("jing_uuid_barcode","w")
    for item in uuid_barcode_dic:
        fout.write(item+"\t"+uuid_barcode_dic[item]+"\n")
    fout.close()

    return uuid_barcode_dic

def uuid_Sample_all():  #the dictionary of all Aliquote key=uuid value=barcode
    localfile = "uuid_Sample.json"
    current=0
    if os.path.exists(localfile):
        fileTime= time.strftime("%Y-%m-%d",time.localtime(os.path.getmtime(localfile)))
        today= datetime.date.today()
        if fileTime ==str(today):
            current=1
    if not current:
        site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/"
        q="json?elementType=Sample"
        os.system("wget -nv -O "+localfile +" "+site+q)
        
    iHandle = open(localfile,"r")
    J = json.loads(iHandle.read())
    iHandle.close()

    dic={}
    for item in J["tcgaElement"]:
        address = item["@href"]
        uuid = string.split(address,"/")[-1]
        dic[uuid]=0

    uuid_barcode_dic={}

    """
    os.system("wget -O uuid_barcode "+uuid_barcode_txt)
    fin = open("uuid_barcode","r")
    for line in fin.readlines():
        uuid, barcode =string.split(line,'\t')[0:2]
        if not uuid_barcode_dic.has_key(uuid):
            uuid_barcode_dic[uuid]=barcode
    fin.close()
    """
    
    if os.path.exists("jing_uuid_barcode_sample"):
        fin = open("jing_uuid_barcode_sample","r")
        for line in fin.readlines():
            line =line[:-1]
            if line =="":
                continue
            uuid, barcode =string.split(line,'\t')[0:2]
            if not uuid_barcode_dic.has_key(uuid):
                uuid_barcode_dic[uuid]=barcode
            uuid_barcode_dic[uuid]=barcode
        fin.close()
        
    for item in dic:
        if item not in uuid_barcode_dic:
            c=1
            while c:
                try:
                    site = "https://tcga-data.nci.nih.gov/uuid/uuidws/metadata/"
                    q="json/uuid/"+item
                    os.system("wget -nv -O tmp "+ site+q)
                    iHandle = open("tmp","r")
                    J = json.loads(iHandle.read())
                    iHandle.close()
                    uuid_barcode_dic[item]=J["tcgaElement"]["barcodes"]["barcode"]
                    break
                except ValueError:
                    c=c +1
                    if c>3:
                        break
    os.system("cp jing_uuid_barcode_sample jing_uuid_barcode_sample_BK")
    fout= open("jing_uuid_barcode_sample","w")
    for item in uuid_barcode_dic:
        fout.write(item+"\t"+uuid_barcode_dic[item]+"\n")
    fout.close()

    return uuid_barcode_dic

UUID_NORMAL_CELLLINE=uuid_normal_cellline()
UUID_CELLLINE=uuid_cellline()
