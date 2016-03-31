import string,sys, os

def getValue(info,key):
    if key in info:
        return info[key]
    else:
        return ""

def parseInfo(info):
    dic ={}
    for item in string.split(info,';'):
        key, value =string.split(item,"=")
        dic[key]=value
    return dic

def ANNheader():
    ANNheader = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' 
    r={}
    items= string.split(ANNheader,'|')
    for i in range(0,len(items)):
        r[string.strip(items[i])]=i
    return r

def parseSnpEffANN (alt, ANN , header):
    for item in string.split(ANN,','):
        if string.strip(item) =="":
            continue
        data = string.split(item,'|')
        if data[header['Allele']] != alt:
            continue
        Annotation =  data[header['Annotation']]
        impact = data[header['Annotation_Impact']]
        gene = data[header['Gene_Name']]
        Feature_Type = data[header['Feature_Type']]
        aa = data[header['HGVS.p']]
        distance = data[header['Distance']]        
        return [Annotation, impact, gene, Feature_Type, aa, distance]

def parseVCF(sample, file, header, fout):
    fin = open(file, 'r')
    for line in fin.readlines():
        if line[0]=="#":
            continue
        data = string.split(line,'\t')
        chr =data[0]
        pos =data[1]
        ref =data[3]
        alt =data[4]
        assert(ref!=".")
        assert(ref!="-")
        assert(alt!=".")
        assert(alt!="-")
        filter = data[6]
        infoLine = data[7]
        if (filter =="LOWSUPPORT"):
            continue
        if chr[:3]!=string.lower("chr"):
            chr = "chr"+chr
        else:
            chr = "chr"+chr[3:]

        start = pos
        end = str(int(pos)+len(ref))
        info = parseInfo(infoLine)
        medianVAF =getValue (info, 'medianVAF')
        ANN =getValue (info, 'ANN')
        Annotation, impact, gene, Feature_Type, aa, distance = parseSnpEffANN(alt, ANN, header)
        if Annotation =="sequence_feature":
            Annotation =""
            gene =""
        if Annotation =="intergenic_region":
            gene =""
        if Annotation =="intron_variant":
            gene =""
        if Annotation == "TF_binding_site_variant":
            Annotation = Feature_Type +" "+Annotation
        outputlist= [sample, chr, start, end, ref, alt, gene, Annotation, medianVAF, aa, impact, distance]
        fout.write(string.join(outputlist, '\t')+'\n')
    fin.close()

def findFiles(inDir, endPattern):
    retfiles =[]
    for root, dirs, files in os.walk(inDir):
        for file in files:
            if file[-len(endPattern):]!=endPattern:
                continue
            print(os.path.join(root, file))
            retfiles.append(os.path.join(root, file))
    return retfiles

if len(sys.argv[:])!= 3:
    print "python parseDirVCF.py inDir outputfile"
    sys.exit()

inDir = sys.argv[1]
if inDir[-1]=="/":
    inDir = inDir[:-1]
files =  findFiles (inDir,  ".eff.vcf")
header = ANNheader()

fout = open(sys.argv[2],'w')
outputList=["sample", "chr", "start","end","reference","alt", "gene", "effect","DNA_VAF", "Amino_Acid_Change","impact","distance"]
fout.write(string.join(outputList,'\t')+'\n')

for file in files:
    sample = string.split(os.path.basename(file),'.')[0]  # pcawg specific
    print sample
    parseVCF(sample, file, header, fout)
fout.close()
