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

uuid_barcode_txt= "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/metadata/metadata.current.txt"

# all sorted by short titles
cancerOfficial={
    "LAML":"acute myeloid leukemia",
    "ACC":"adrenocortical carcinoma",
    "BLCA":"bladder urothelial carcinoma",
    "LGG":"brain lower grade glioma",
    "BRCA":"breast invasive carcinoma",
    "CESC":"cervical squamous cell carcinoma and endocervical adenocarcinoma",
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
    "LUNG":"lung cancer",
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
    "CHOL":"cholangiocarcinoma",
    "LCLL":"chronic lymphocytic leukemia",
    "LCML":"chronic myelogenous leukemia",
    "MESO":"mesothelioma",
    "MISC":"miscellaneous",
    "PCPG":"pheochromocytoma and paraganglioma",
    "UCS":"uterine carcinosarcoma",
    "UVM":"uveal melanoma",
    "TGCT":"testicular germ cell tumor",
    "PANCAN":"pan-cancer"
    }

cancerHumanReadable={
    "LAML":"Acute Myeloid Leukemia",
    "ACC":"Adrenocortical Cancer",
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
    "LUNG":"Lung Cancer",
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
    "CHOL":"Bile Duct Cancer",
    "LCLL":"Chronic Lymphocytic Leukemia",
    "LCML":"Chronic Myeloid Leukemia",
    "MESO":"Mesothelioma",
    "MISC":"Miscellaneous",
    "PCPG":"Pheochromocytoma & Paraganglioma",
    "UCS":"Uterine Carcinosarcoma",
    "UVM":"Ocular melanomas",
    "TGCT":"Testicular Cancaer",
    "PANCAN":"Pan-Cancer"
    }

anatomical_origin ={
    "LAML":["White blood cell"],
    "ACC":["Adrenal gland"],
    "BLCA":["Bladder"],
    "LGG":["Brain"],
    "BRCA":["Breast"],
    "CESC":["Cervix"],
    "CNTL":[],
    "COAD":["Colon"],
    "COADREAD":["Colon and Rectum"],
    "ESCA":["Esophagus"],
    "GBM":["Brain"],
    "HNSC":["Head and Neck region"],
    "KICH":["Kidney"],
    "KIRC":["Kidney"],
    "KIRP":["Kidney"],
    "LIHC":["Liver"],
    "LUAD":["Lung"],
    "LUNG":["Lung"],
    "LUSC":["Lung"],
    "DLBC":["Lymphatic tissue"],
    "LNNH":["Lymphatic tissue"],
    "OV":["Ovary"],
    "PAAD":["Pancreas"],
    "PRAD":["Prostate"],
    "READ":["Rectum"],
    "SARC":["Soft tissue", "Bone"],
    "SKCM":["Skin"],
    "STAD":["Stomach"],
    "THCA":["Thyroid Gland"],
    "UCEC":["Uterus"],
    "CHOL":["Bile duct"],
    "LCLL":["Blood"],
    "LCML":["Blood"],
    "MESO":["Lining of body cavities"],
    "MISC":[],
    "PCPG":["Paraganglia"],
    "UCS":["Uterus"],
    "UVM":["Eye"],
    "TGCT":["Testis"],
    "PANCAN":[]
    }

tags ={
    "LAML":[],
    "ACC":[],
    "BLCA":[],
    "LGG":["nervous system"],
    "BRCA":[],
    "CESC":[],
    "CNTL":[],
    "COAD":[],
    "COADREAD":[],
    "ESCA":[],
    "GBM":["nervous system"],
    "HNSC":[],
    "KICH":[],
    "KIRC":[],
    "KIRP":[],
    "LIHC":[],
    "LUAD":["non-small cell lung cancer"],
    "LUNG":["non-small cell lung cancer"],
    "LUSC":["non-small cell lung cancer"],
    "DLBC":[],
    "LNNH":[],
    "OV":[],
    "PAAD":[],
    "PRAD":[],
    "READ":[],
    "SARC":[],
    "SKCM":[],
    "STAD":["gastric cancer"],
    "THCA":[],
    "UCEC":[],
    "CHOL":[],
    "LCLL":[],
    "LCML":[],
    "MESO":[],
    "MISC":[],
    "PCPG":[],
    "UCS":[],
    "UVM":[],
    "TGCT":[],
    "PANCAN":["non-small cell lung cancer","gastric cancer"]
    }

cancerGroupTitle={
    "LAML":"acute myeloid leukemia",
    "BLCA":"bladder urothelial carcinoma",
    "LGG":"brain lower grade glioma",
    "BRCA":"breast invasive carcinoma",
    "CESC":"cervical & endocervical cancer",
    "CNTL":"controls",
    "COAD":"colon adenocarcinoma",
    "COADREAD":"colon & rectum adenocarcinoma",
    "ESCA":"esophageal carcinoma",
    "GBM":"glioblastoma multiforme",
    "HNSC":"head & neck squamous cell carcinoma",
    "KICH":"kidney chromophobe",
    "KIRC":"kidney clear cell carcinoma",
    "KIRP":"kidney papillary cell carcinoma",
    "LIHC":"liver hepatocellular carcinoma",
    "LUAD":"lung adenocarcinoma",
    "LUNG":"lung cancer",
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
    "UCS":"uterine carcinosarcoma",
    "CHOL":"cholangiocarcinoma",
    "UVM":"uveal Melanoma",
    "TGCT":"testicular germ cell tumor",
    "PANCAN":" Pan-Cancer"
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
    "sample_type_id":"category",
    "_OS_IND":"category",
    "_RFS_IND":"category",
    "_EVENT":"category",
    "":""
    }

featureStateOrder={
    "er_level_cell_percentage_category":{"BRCA":
                                         ["<10%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-99%"]},
    "breast_carcinoma_estrogen_receptor_status":{"BRCA":
                                                 ["Negative", "Indeterminate",  "Positive"]},
    "breast_carcinoma_progesterone_receptor_status":{"BRCA":
                                                 ["Negative", "Indeterminate",  "Positive"]},
    "her2_erbb_pos_finding_cell_percent_category":{"BRCA":
                                                   ["<10%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-99%"]},
    "distant_metastasis_pathologic_spread":{"ALL":
                                            ['MX', 'M0', 'M1', 'M1a', 'M1b']},
    "histological_type":{"LGG":
                         ["Oligodendroglioma","Oligoastrocytoma","Astrocytoma"]},
    "sample_type":{"ALL":
                   ["Primary Tumor",
                    "Primary solid Tumor",
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
    "tobacco_smoking_history":{"ALL":["Lifelong Non-smoker",
                                                "Current reformed smoker for > 15 years",
                                                "Current reformed smoker for < or = 15 years",
                                                "Current smoker"]
                                         },
    "vital_status":{"ALL": ["LIVING","DECEASED"]},
    "mononucleotide_and_dinucleotide_marker_panel_analysis_status": {"ALL": [ "MSI-H", "MSI-L", "Indeterminate","MSS"]},
    "neoplasm_histologic_grade": {"ALL" :["Low Grade","High Grade" ]},
    "":""
    }

featurePriority={
    "BLCA": {"sample_type":"1",
             "_PANCAN_Cluster_Cluster_PANCAN":"2",
             "neoplasm_histologic_grade":"3",
             "gender":"4"
             },
    "BRCA": {"sample_type":"1",
             "PAM50 mRNA_nature2012":"2",
             "ER Status_nature2012":"3",
             "PR Status_nature2012":"4",
             "HER2 Final Status_nature2012":"5"
             },
    "COAD": {"sample_type":"1",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"2",
             "histological_type":"3",
             "CIMP":"4"
             },
    "COADREAD":{"sample_type":"1",
                "_PANCAN_DNAMethyl_PANCAN":"2",
                "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"3",
                "cohort":"4"
                },
    "CESC": {"sample_type":"1",
             "neoplasm_histologic_grade":"2",
             "histological_type":"3",
             "_TIME_TO_EVENT":"4"
             },
    "ESCA": {"sample_type":"1"
             },
    "DLBC": {"sample_type":"1"
             },
    "GBM":  {"sample_type":"1",
             "GeneExp_Subtype":"2",
             "G_CIMP_STATUS":"3",
             "_PANCAN_DNAMethyl_PANCAN":"4"
             },
    "HNSC": {"sample_type":"1",
             "hpv_status_by_p16_testing":"2",
             "hpv_status_by_ish_testing":"3",
             "gender":"4"
             },
    "KICH": {"sample_type":"1",
             "person_neoplasm_cancer_status":"2"
             },
    "KIRC": {"sample_type":"1",
             "neoplasm_histologic_grade":"2",
             "person_neoplasm_cancer_status":"3"
             },
    "KIRP": {"sample_type":"1",
             "person_neoplasm_cancer_status":"2"
             },
    "LIHC": {"sample_type":"1",
             "neoplasm_histologic_grade":"2",
             },
    "LAML": {"sample_type":"1",
             "acute_myeloid_leukemia_calgb_cytogenetics_risk_category":"2",
             "_PANCAN_DNAMethyl_LAML":"3",
             "_PANCAN_mirna_LAML":"4",
             "cytogenetic_abnormality":"5"
             },
    "LGG":  {"sample_type":"1",
             "histological_type":"2",
             "neoplasm_histologic_grade":"3",
             "_TIME_TO_EVENT":"4"
             },
    "LUAD": {"sample_type":"1",
             "Expression_Subtype":"2",
             "EGFR":"3",
             "KRAS":"4",
             "tobacco_smoking_history":"5"
             },
    "LUNG": {"sample_type":"1",
             "cohort":"2",
             "_PANCAN_UNC_RNAseq_PANCAN_K16":"3",
             "_PANCAN_DNAMethyl_PANCAN":"4"
             },
    "LUSC": {"sample_type":"1",
             "_PANCAN_UNC_RNAseq_PANCAN_K16":"2",
             "_PANCAN_DNAMethyl_PANCAN":"3",
             "tobacco_smoking_history":"4"
             },
    "KICH": {"sample_type":"1"
             },
    "OV": {"sample_type":"1",
           "neoplasm_histologic_grade":"2"
           },
    "PAAD": {"sample_type":"1"
           },
    "PRAD": {"sample_type":"1",
             "gleason_score":"2",
             "residual_tumor":"3",
             "psa_value":"4"
             },
    "READ": {"sample_type":"1",
             "histological_type":"2",
             "_PANCAN_DNAMethyl_PANCAN":"3",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"4"
             },
    "SARC": {"sample_type":"1",
             "histological_type":"2",
             "tumor_tissue_site":"3",
             "gender":"4"
             },
    "SKCM": {"sample_type":"1",
             "postoperative_rx_tx":"2",
             "melanoma_clark_level_value":"3",
             "gender":"4"
             },
    "STAD": {"sample_type":"1",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"2",
             "neoplasm_histologic_grade":"3",
             "h_pylori_infection":"4"
             },
    "THCA": {"sample_type":"1",
             "pathologic_stage":"2",
             "histological_type":"3",
             "gender":"4"
             },
    "UCEC": {"sample_type":"1",
             "histological_type":"2",
             "mononucleotide_and_dinucleotide_marker_panel_analysis_status":"3",
             "neoplasm_histologic_grade":"4"
             },
    "ACC":{"sample_type":"1"
             },
    "LCLL":{"sample_type":"1"
            },
    "LCML":{"sample_type":"1"
             },
    "MESO":{"sample_type":"1"
            },
    "MISC":{"sample_type":"1"
            },
    "PCPG":{"sample_type":"1"
            },
    "TGCT":{"sample_type":"1"
            },
    "UCS":{"sample_type":"1"
           },
    "UCEC": {"sample_type":"1",
             "cohort":"2",
             "_PANCAN_UNC_RNAseq_PANCAN_K16":"3",
             "_PANCAN_DNAMethyl_PANCAN":"4"
             },
    "PANCAN": {"sample_type":"1",
               "cohort":"2",
               "_PANCAN_UNC_RNAseq_PANCAN_K16":"3",
               "gender":"4"
             }
    }

clinDataDesc ="The accompanied clinical data are downloaded from TCGA data coordination center. They are the patient, sample, follow up, biospecimen tumor sample, and auxiliary clinical information. Names of the clinical fields were \"curated\" by UCSC to be more readable."

clinDataDescBRCA="There are now clinical features from the Supplemental table 1 from the Nature 2012 paper (pubmed:23000897) including HER2, PR and ER status as well as PAM50 calls based off the Agilent 244K custom gene expression microarrays. We have also included preliminary PAM50 calls based on the Illumina HiSeq RNA Sequencing platform from the TCGA Analysis Working Group (AWG). Note that these calls are not final and are subject to change."

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
