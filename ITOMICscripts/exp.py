#https://github.com/jamescasbon/PyVCF
import sys, os, string
import vcf

vcffile = "DNARNA"
sample = "TUMOR_RNA"
outdir ="./"

def findSamples(vcffile):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    return vcf_reader.samples

def exp_data (vcffile, sample, outdir, sampleLabel):
    parseIsoformGeneLabel = lambda x: string.split(x,'/')
    EXP_I_PCT_dic ={}
    EXP_I_CNT_dic ={}
    EXP_I_TPM_dic ={}
    EXP_I_FPKM_dic ={}
    EXP_G_CNT_dic ={}
    EXP_G_TPM_dic ={}
    EXP_G_FPKM_dic ={}

    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    for record in vcf_reader:
        if record.ALT[0].type == "EXP":
            data = parseIsoformGeneLabel(record.INFO["ISO"])
            isoform = data[0]
            hugos= data[1:]
            #if record.INFO.has_key("CANONICAL"):
            #    print record.INFO["CANONICAL"]

            EXP_I_PCT = record.genotype(sample)["EXP_I_PCT"]
            EXP_I_CNT = record.genotype(sample)["EXP_I_CNT"]
            EXP_I_TPM = record.genotype(sample)["EXP_I_TPM"]
            EXP_I_FPKM = record.genotype(sample)["EXP_I_FPKM"]
            EXP_G_CNT = record.genotype(sample)["EXP_G_CNT"]
            EXP_G_TPM = record.genotype(sample)["EXP_G_TPM"]
            EXP_G_FPKM = record.genotype(sample)["EXP_G_FPKM"]

            EXP_I_PCT_dic[isoform]=EXP_I_PCT
            EXP_I_CNT_dic[isoform]=EXP_I_CNT
            EXP_I_TPM_dic[isoform]=EXP_I_TPM
            EXP_I_FPKM_dic[isoform]=EXP_I_FPKM
            for hugo in hugos:
                EXP_G_CNT_dic[hugo]=EXP_G_CNT
                EXP_G_TPM_dic[hugo]=EXP_G_TPM
                EXP_G_FPKM_dic[hugo]=EXP_G_FPKM
    fin.close()

    output_dic (sampleLabel, "IsoPct", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_I_PCT")
    output_dic (sampleLabel, "expected_count", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_I_CNT")
    output_dic (sampleLabel, "tpm", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_I_TPM")
    output_dic (sampleLabel, "fpkm", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_I_FPKM")
    output_dic (sampleLabel, "expected_count", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_G_CNT")
    output_dic (sampleLabel, "tpm", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_G_TPM")
    output_dic (sampleLabel, "fpkm", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_G_FPKM")

def output_dic (sample, unit, dic, file):
    fout=open(file,'w')
    fout.write(sample +"\t"+unit+"\n")
    for key in dic.keys():
        fout.write(key +"\t"+str(dic[key])+"\n")
    fout.close()

samples =  findSamples(vcffile)
if sample not in samples:
    print "bad sample name"
    sys.exi()

exp_data(vcffile, sample, outdir)
