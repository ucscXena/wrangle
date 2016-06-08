#https://github.com/jamescasbon/PyVCF

import string, vcf

def checkSample(vcffile, sample):
    samples =  findSamples(vcffile)
    if sample not in samples :
        return False
    return True

def findSamples(vcffile):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    return vcf_reader.samples

def findAssembly (vcffile):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    return vcf_reader.metadata["reference"]

def findFileMetaData (vcffile):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    return vcf_reader.metadata

def findSampleMetaData (vcffile, sample):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    for item in vcf_reader.metadata['SAMPLE']:
        if item["ID"] == sample:
            return item
    return None

def printRecordInfo(vcf_record, vcf_reader):#vcf_record
    infos = vcf_reader.infos
    for key in vcf_record.INFO.keys():
        if key in infos:
            print key, vcf_record.INFO[key], infos[key]
        else:
            print key, vcf_record.INFO[key]

#http://www.cravat.us/help.jsp
CRAVAT_SO_code ={
    "SY":"synonymous_variant",
    "SL":"stop_lost",
    "SG":"stop_gained",
    "MS":"missense_variant",
    "II":"inframe_insertion",
    "FI":"Frame_Shift_Ins",
    "ID":"inframe_deletion",
    "FD":"Frame_Shift_Del",
    "CS":"Complex Substitution"
}
