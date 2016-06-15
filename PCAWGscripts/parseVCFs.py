#https://github.com/jamescasbon/PyVCF
import sys, os, string, copy
import vcf, json, uuid
import urllib2
sys.path.insert(0, os.path.dirname(sys.argv[0])+"/../support/")
import xenaVCF
import probMap_genePred

sampleTumor = "TUMOUR"
sampleNormal = "NORMAL"
start_sv_padding = 2000
promoter_padding = 1000
end_padding = 0
ann_url = "https://reference.xenahubs.net/download/gencode_good_hg19"

#annotation a SV : if SV cut within a gene +- Nbp
def annotate_SV(chr, start, end, start_padding, end_padding, annDic):
    genes = []
    #gene annotation
    for hugo in annDic[chr].keys():
        for item in annDic[chr][hugo]:
            chrHugo = item['chr']
            startHugo = item ['start']
            endHugo = item['end']
            strand = item['strand']
            if  chrHugo != chr:
                continue
            if strand =="+":
                if start <= endHugo  + end_padding  and end  > startHugo - start_padding:
                    if hugo not in genes:
                        genes.append(hugo)
            if strand =="-":
                if start <= endHugo +start_padding  and end  > startHugo - end_padding:
                    if hugo not in genes:
                        genes.append(hugo)
    return genes

def annotate_SNV_extended_exon (chr, start, end, start_padding, end_padding, annDic):
    genes = []
    #gene annotation
    for hugo in annDic[chr].keys():
        if hugo in genes:
            continue
        for item in annDic[chr][hugo]:
            chrHugo = item['chr']
            if  chrHugo != chr:
                continue

            startHugo = item ['start']
            endHugo = item['end']
            strand = item['strand']

            if strand =="+":
                if start <= startHugo - start_padding  or end  > endHugo + end_padding : # outside the gene region
                    continue

            if strand =="-":
                if start <= startHugo - end_padding  or end  > endHugo + start_padding : # outside the gene region
                    continue

            exonStarts = copy.deepcopy(item['exonStarts'])
            exonEnds = copy.deepcopy(item['exonEnds'])
            exonCount = item['exonCount']

            if strand =="+":
                exonStarts[0]=exonStarts[0] - start_padding
                exonEnds[-1]=exonEnds[-1] + end_padding

            if strand =="-":
                exonStarts[0] = exonStarts[0] - end_padding
                exonEnds[-1] = exonEnds[-1] + start_padding

            for i in range (0, exonCount):
                if start <= exonEnds[i]  and end  > exonStarts[i]:
                    if hugo not in genes:
                        genes.append(hugo)
                        break

            if hugo in genes:
                break

    return genes

def parse_BND (vcffile, start_padding, end_padding, annDic):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    ret_data =[]
    for record in vcf_reader:
        #print record.ALT, len(record.ALT), record.ALT[0].type 
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
            genes =  (chr, start, end, start_padding, end_padding, annDic)
            for gene in genes:
                data={}
                data['chr']= chr
                data['start']= start
                data['end']= end
                data['ref']= ref
                data['alt']= alt
                data['gene']=gene
                data['effect']="SV"
                ret_data.append(data)
    fin.close()
    return ret_data


def cavat_annotate(chr, start, ref, alt, count):
    coding_gene = None
    effect = ""
    aa = ""

    #if count >10 :
    #    return count, coding_gene, effect, aa
    try:
        url = "http://www.cravat.us/rest/service/query?mutation="+chr+"_"+str(start)+"_+_"+ref+ "_" + str(alt)
        cravat = urllib2.urlopen(url)
        ann=  json.loads(cravat.read())
        cravat.close()

        if ann["HUGO symbol"] and ann["HUGO symbol"]!="Non-Coding":
            coding_gene = ann["HUGO symbol"]
                
        #CRAVAT analysis
        #http://www.cravat.us/help.jsp
        if ann["Sequence ontology"]!="":
            effect = ann["Sequence ontology"]
            if effect in xenaVCF.CRAVAT_SO_code:
                effect = xenaVCF.CRAVAT_SO_code[effect]

        if ann["Sequence ontology protein change"]!="":
            aa = ann["Sequence ontology protein change"]

        return count, coding_gene, effect, aa

    except:
        return cavat_annotate(chr, start, ref, alt, count +1)


def parse_SNV (vcffile, start_padding, end_padding, annDic):
    fin=open(vcffile, 'r')
    vcf_reader = vcf.Reader(fin)
    ret_data =[]
    for record in vcf_reader:
        #print record.ALT, len(record.ALT), record.ALT[0].type 
        type =  record.ALT[0].type
        if type in ["SNV"]:
            # assuming no filter is passing
            if record.FILTER and len(record.FILTER)!=0: 
                #print record, record.FILTER
                continue

            chr = "chr" + record.CHROM
            start = record.POS
            end = start + len(record.REF) -1
            ref = record.REF
            alt = record.ALT[0]
            effect =""
            aa =""
            try:
                VAF = record.INFO["VAF"]
            except KeyError:
                VAF=""

            #promoter, UTR, noncoding 
            all_genes = annotate_SNV_extended_exon  (chr, start, end, start_padding, end_padding, annDic)
            coding_genes =[]

            if len(all_genes) ==0:
                continue

            count =0
            r = cavat_annotate(chr, start, ref, alt, count)
            count, coding_gene, effect, aa = r

            if coding_gene != '' :
                coding_genes = [coding_gene]
            if count >1:
                url = "http://www.cravat.us/rest/service/query?mutation="+chr+"_"+str(start)+"_+_"+ref+ "_" + str(alt)
                print count, url

            for gene in coding_genes:
                data={}
                data['chr']= chr
                data['start']= start
                data['end']= end
                data['ref']= ref
                data['alt']= alt
                data['gene']=gene
                data['aa']= aa
                data['effect']=effect
                data['VAF']=VAF
                ret_data.append(data)

            for gene in all_genes:
                if gene not in coding_genes:
                    data={}
                    data['chr']= chr
                    data['start']= start
                    data['end']= end
                    data['ref']= ref
                    data['alt']= alt
                    data['gene']=gene
                    data['aa']= ""
                    data['effect']="unknown"
                    data['VAF']=VAF
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
        fout.write(sample)
        fout.write('\t'+ item['chr'])
        fout.write('\t'+ str(item['start']))
        fout.write('\t'+ str(item['end']))
        fout.write('\t'+ str(item['ref']))
        fout.write('\t'+ str(item['alt']))
        fout.write('\t'+ str(item['gene']))
        fout.write('\t'+ item['effect'])
        fout.write('\t'+ str(item['VAF']))
        fout.write('\t'+str(item['aa']))
        fout.write('\n')

def cleanSVPCAWGvcf(file): #stupid
    output = str(uuid.uuid4())
    os.system("grep -v ^##contig "+ file +" |grep -v ^##pcawg > "+output)
    return output

if len(sys.argv[:])!=4:
    print "python parseVCFs.py listFile(one_filename_perline) dataType(BND,SNV) outputfile"
    sys.exit()

flist = open(sys.argv[1],'r')
dataType = sys.argv[2]
fout = open(sys.argv[3],'w')

stream = urllib2.urlopen(ann_url)
annDic = probMap_genePred.parseGenePredToGene(stream)
stream.close()

for infile in flist.readlines():
    infile = infile[:-1]
    if dataType =="BND":
        vcffile = cleanSVPCAWGvcf(infile)
        if not xenaVCF.checkSample(vcffile, sampleTumor):
            print sampleTumor, "bad sample name" 
            os.system("rm -f "+ vcffile)
            continue

        tumorMetaData = xenaVCF.findSampleMetaData(vcffile,sampleTumor)
        sampleLabel = tumorMetaData['SampleName']
        if dataType =="BND":
            xenaRecords=  parse_BND (vcffile, start_sv_padding, end_padding,  annDic)
            outputputMutationVector (sampleLabel, xenaRecords, fout)

    elif dataType == "SNV":
        sampleLabel = string.split(os.path.basename(infile),'.')[0] # stupid PCAWG SNV VCFS
        print sampleLabel
        xenaRecords=  parse_SNV (infile, promoter_padding, end_padding, annDic)
        outputputMutationVector (sampleLabel, xenaRecords, fout)

flist.close()
fout.close()

#assembly = xenaVCF.findAssembly (vcffile)
#print assembly

#ploidy = tumorMetaData['Ploidy']
#print ploidy




