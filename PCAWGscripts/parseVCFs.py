#https://github.com/jamescasbon/PyVCF
import sys, os, string
import vcf, json, uuid
import urllib2

sys.path.insert(0, os.path.dirname(sys.argv[0])+"../support/")

import xenaVCF
import probMap_genePred

#annotation a SV : if SV cut within a gene +- Nbp
def annotate_SV(chr, start, end, annDic):
    genes = []
    #gene annotation
    for hugo in annDic[chr].keys():
        for item in annDic[chr][hugo]:
            chrHugo = item['chr']
            startHugo = item ['start']
            endHugo = item['end']
            if  chrHugo != chr:
                continue
            if start <= endHugo and end > startHugo:
                if hugo not in genes:
                    genes.append(hugo)
    return genes

def parse_BND (vcffile, annDic):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    ret_data =[]
    for record in vcf_reader:
        #print record.ALT, len(record.ALT), record.ALT[0].type 
        data={}
        type =  record.ALT[0].type
        if type in ["BND"]:
            # assuming no filter is passing
            if len(record.FILTER)!=0: 
                print record, record.FILTER
                
            #print record.genotype(sampleTumor)
            #print record.genotype(sampleNormal)
            #printRecordInfo(record, vcf_reader) 
                   
            chr = "chr" + record.CHROM
            start = record.POS
            end = start + len(record.REF) -1
            ref = record.REF
            alt = record.ALT[0]
            genes = annotate_SV(chr, start, end, annDic)
            for gene in genes:
                data['chr']= chr
                data['start']= start
                data['end']= end
                data['ref']= ref
                data['alt']= alt
                data['gene']=gene
                ret_data.append(data)
    fin.close()
    return ret_data

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

def outputputMutationVector (sample, dataList, fout):
    for item in dataList:
        fout.write(sample+'\t')
        fout.write(item['chr']+'\t')
        fout.write(str(item['start'])+'\t')
        fout.write(str(item['end'])+'\t')
        fout.write(str(item['ref'])+'\t')
        fout.write(str(item['alt'])+'\t')
        fout.write(str(item['gene'])+'\t')
        fout.write('\n')

def cleanPCAWGvcf(file): #stupid
    output = str(uuid.uuid4())
    os.system("grep -v ^##contig "+ file +" |grep -v ^##pcawg > "+output)
    return output


if len(sys.argv[:])!=4:
    print "python parseVCFs.py listFile(one_filename_perline) dataType(BND) outputfile"
    sys.exit()

flist = open(sys.argv[1],'r')
dataType = sys.argv[2]
fout = open(sys.argv[3],'w')

sampleTumor = "TUMOUR"
sampleNormal = "NORMAL"
ann_url = "https://reference.xenahubs.net/download/gencode_good_hg19"
stream = urllib2.urlopen(ann_url)
annDic = probMap_genePred.parseProbeMapToGene(stream)
stream.close()


for infile in flist.readlines():
    infile = infile[:-1]
    vcffile = cleanPCAWGvcf(infile)
    if not xenaVCF.checkSample(vcffile, sampleTumor):
        print sampleTumor, "bad sample name" 
        os.system("rm -f "+ vcffile)
        continue

    tumorMetaData = xenaVCF.findSampleMetaData(vcffile,sampleTumor)
    sampleLabel = tumorMetaData['SampleName']
    if dataType =="BND":
        xenaRecords=  parse_BND (vcffile, annDic)
        outputputMutationVector (sampleLabel, xenaRecords, fout)
        os.system("rm -f "+ vcffile)

flist.close()
fout.close()

#assembly = xenaVCF.findAssembly (vcffile)
#print assembly

#ploidy = tumorMetaData['Ploidy']
#print ploidy




