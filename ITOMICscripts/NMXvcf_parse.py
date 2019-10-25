#https://github.com/jamescasbon/PyVCF
import sys, os, string
import vcf, json

refGeneTranscriptParse = lambda x: string.split(x,".")[0]

def findSamples(vcffile):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    fin.close()
    return vcf_reader.samples

def EXP_data (vcffile, sample, outdir, sampleLabel):
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
        type = record.ALT[0].type
        if type == "EXP":
            data = parseIsoformGeneLabel(record.INFO["ISO"])
            isoform = refGeneTranscriptParse(data[0])
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

    output_dic (sampleLabel, "IsoPct", "transcript expression RNAseq", EXP_I_PCT_dic, outdir+ sampleLabel+"_EXP_I_PCT")
    output_dic (sampleLabel, "expected_count",  "transcript expression RNAseq", EXP_I_CNT_dic, outdir+ sampleLabel+"_EXP_I_CNT")
    output_dic (sampleLabel, "tpm",  "transcript expression RNAseq", EXP_I_TPM_dic, outdir+ sampleLabel+"_EXP_I_TPM")
    output_dic (sampleLabel, "fpkm",  "transcript expression RNAseq", EXP_I_FPKM_dic, outdir+ sampleLabel+"_EXP_I_FPKM")
    output_dic (sampleLabel, "expected_count",  "gene expression RNAseq", EXP_G_CNT_dic, outdir+ sampleLabel+"_EXP_G_CNT")
    output_dic (sampleLabel, "tpm",  "gene expression RNAseq", EXP_G_TPM_dic, outdir+ sampleLabel+"_EXP_G_TPM")
    output_dic (sampleLabel, "fpkm",  "gene expression RNAseq", EXP_G_FPKM_dic, outdir+ sampleLabel+"_EXP_G_FPKM")


def CNV_data (vcffile, sampleTumor, sampleNormal, outdir, sampleLabel):
    output = outdir+ sampleLabel+"_CNV"
    fout = open(output,'w')
    fout.write(string.join(["chr", "start", "end", "TUMOR_CNV_RC", "TUMOR_AF", "TUMOR_AFN"],'\t')+'\n')

    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    for record in vcf_reader:
        #print record.ALT, len(record.ALT), record.ALT[0].type 

        type =  record.ALT[0].type
        if type in ["CNV"]:
            if len(record.FILTER)!=0:
                print (record, record.FILTER)

            #print record.genotype(sampleTumor)
            #print record.genotype(sampleNormal)

            chr = record.CHROM
            start = record.POS
            size = record.INFO['SVLEN']
            end = record.POS + size -1
            TUMOR_CNV_RC = record.genotype(sampleTumor)["CNV_RC"]
            TUMOR_AF= record.genotype(sampleTumor)["CNV_AF"]
            TUMOR_AFN = record.genotype(sampleTumor)["CNV_AFN"]
            
            if TUMOR_AFN==0:
                TUMOR_AF=""
                TUMOR_AFN=""
            fout.write( string.join([chr, str(start), str(end), str(TUMOR_CNV_RC), str(TUMOR_AF), str(TUMOR_AFN)],'\t')+'\n')

    fin.close()


def output_dic (sample, unit, dataSubType, dataDic, file):
    fout=open(file,'w')
    fout.write( "id\t"+sample+"\n")
    for key in dataDic.keys():
        fout.write(key +"\t"+str(dataDic[key])+"\n")
    fout.close()

    fout=open(file+".json",'w')
    j={}
    j["unit"]=unit
    j["type"]="genomicMatrix"
    j["dataSubType"]= dataSubType
    fout.write( json.dumps( j, indent=-1 ) )
    fout.close()

if len(sys.argv[:])!=4:
    print ("python NMXvcf_parse.py vcf dataType(EXP,CNV) sampleLabel")
    sys.exit()

vcffile = sys.argv[1]
dataType = sys.argv[2]
sampleLabel = sys.argv[3]
outdir ="./"

samples =  findSamples(vcffile)
print (samples)

sampleTumorRNA = "TUMOR_RNA"
sampleTumor = "TUMOR"
sampleNormal = "NORMAL"

if sampleTumor not in samples or sampleNormal not in samples:
    print ("bad sample name")
    sys.exit()

if sampleTumorRNA not in samples:
    print ("bad sample name")
    sys.exit()

if dataType == "EXP":
    EXP_data(vcffile, sampleTumorRNA, outdir, sampleLabel)
elif dataType =="CNV":
    CNV_data(vcffile, sampleTumor, sampleNormal, outdir, sampleLabel)
