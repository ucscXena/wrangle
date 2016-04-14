import sys,string,os
import json,datetime
import math
import copy
import uuid

def RSEM_IsoPct (inDir, outfile):
    PATHPATTERN= ".rsem_isoforms.results"
    valuePOS=7
    LOG2=0
    UNIT = "IsoPct"
    theta=0
    dataSubType="transcript expression RNAseq"
    geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT, dataSubType)
    return

def RSEM_tpm (inDir, outfile):
    PATHPATTERN= ".rsem_isoforms.results"
    valuePOS=4
    LOG2=1
    UNIT = "tpm"
    theta=0.001
    dataSubType="transcript expression RNAseq"
    geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT,dataSubType)
    return

def RSEM_Hugo_TOIL_norm_counts (inDir, outfile):
    PATHPATTERN= ".rsem.genes.norm_counts.hugo.tab"
    valuePOS=1
    LOG2=1
    UNIT = "norm_counts"
    theta=1
    dataSubType="gene expression RNAseq"
    geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT, dataSubType)
    return

def Kallisto_est_counts (inDir, outfile):
    PATHPATTERN= ".abundance.tsv"
    valuePOS=3
    LOG2=1
    theta=1
    UNIT="est_counts"
    dataSubType="transcript expression RNAseq"
    geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT, dataSubType)
    return

def Kallisto_tpm (inDir, outfile):
    PATHPATTERN= ".abundance.tsv"
    valuePOS=4
    LOG2=1
    theta=0.001
    UNIT = "tpm"
    dataSubType="transcript expression RNAseq"
    geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT, dataSubType)
    return

def geneRPKM (inDir, outfile, PATHPATTERN, valuePOS, LOG2, theta, UNIT, dataSubType):
    allSamples={}
    c=0
    dataMatrix=[]
    tmpSamples={}
    genes={}
    oldgenes={}
    files=[]
    ERROR=1
    tmpDir =str(uuid.uuid4())
    #data specific
    RANK=0

    os.system("mkdir "+ tmpDir)

    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)==-1:
            continue
    
        infile = inDir+"/"+file
        #stupid TOIL naming
        sample = string.split(file,".")[0]

        if sample in allSamples:
            print len(allSamples)
            print "ERROR duplicated sample = "+ sample
            continue
                
        p=len(allSamples)
        allSamples[sample]=p
                    

        p=len(tmpSamples)
        tmpSamples[sample]=p
                
        c=c+1
        print c, sample

        process(dataMatrix,tmpSamples,sample,genes, infile,valuePOS,LOG2,theta,250)
        if (c % 250)==0:
            tmpout= tmpDir+"/"+"tmp_"+ str(int(c/250.0))
            ERROR =outputMatrix(dataMatrix, tmpSamples, genes, oldgenes, tmpout)
            assert(not ERROR)

            dataMatrix=[]
            tmpSamples={}
            oldgenes=copy.deepcopy(genes)
            genes ={}
            files.append(tmpout)
                    

    if (c % 250)!=0:
        tmpout=  tmpDir+"/"+"tmp_"+ str(int(c/250.0)+1)
        files.append(tmpout)
        ERROR= outputMatrix(dataMatrix, tmpSamples, genes, oldgenes, tmpout)
        assert(not ERROR)

    #paste all together
    os.system("paste -d \'\' "+string.join(files," ")+" > "+ outfile)

    oHandle = open(outfile+".json","w")    
    J={}
    J["redistribution"]= True
    J["dataProducer"]= "UCSC TOIL RNA-seq recompute"
    J["colNormalization"]=True
    J["type"]= "genomicMatrix" 
    J["version"]= datetime.date.today().isoformat()
    J["wrangler"]= "Xena scripts processed on "+ datetime.date.today().isoformat()
    J["dataSubType"]= dataSubType 
    if LOG2:
        J["unit"]="log2("+UNIT+"+"+str(theta)+")"
        J["wrangling_procedure"]= "Data (file names: *"+PATHPATTERN+") are downloaded, "+ UNIT+" values are extracted, log2(x+"+str(theta)+") transformed, and combined."
    else:
        J["unit"]=UNIT
        J["wrangling_procedure"]= "Data (file names: *"+PATHPATTERN+") are downloaded, "+ UNIT + "values are extracted and combined."

    J["label"] = ""
    J["cohort"] =""

    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()
    
    return

def process(dataMatrix,samples, sample,genes, infile,valuePOS, LOG2, theta, maxLength):
    # one sample a file  
    fin=open(infile,'U')    
    fin.readline()
    for line in fin.readlines():
        data =string.split(line[:-1],"\t")
        hugo = data[0]
        value= data[valuePOS]
        hugo = string.split(hugo,"|")[0]

        if hugo not in genes:
            p=len(genes)
            genes[hugo]=p
            l=[]
            for j in range (0,maxLength):
                l.append("")    
            dataMatrix.append(l)
        
        if value not in ["","null","NULL","Null","NA"]:
            if LOG2:
                value = float(value)
                if value<0:
                    value = ""
                else:
                    value = math.log(float(value+ theta),2)

            x=genes[hugo]
            y=samples[sample]
            dataMatrix[x][y]=value
            
    fin.close()
    return 

def outputMatrix(dataMatrix, samples, genes, oldgenes,outfile):
    #compare genes and oldgenes:
    if oldgenes!={}:
        if len(genes)!=len(oldgenes):
            print "ERROR genes total length is different"
            return 1
        for gene in genes:
            if genes[gene]!=oldgenes[gene]:
                print "ERROR gene order is different",gene
                return 1
            
    fout = open(outfile,"w")
    if oldgenes=={}:
        fout.write("sample")
    for sample in samples:
        fout.write("\t"+sample)
    fout.write("\n")

    for gene in genes:
        if oldgenes=={}:
            fout.write(gene)
        for sample in samples:
            value = dataMatrix[genes[gene]][samples[sample]]
            if value !="":
                value = "%.4f" % (float(value))
                fout.write("\t"+value)
            else:
                fout.write("\tNA")
        fout.write("\n")
    fout.close()
    return 0

if len(sys.argv[:])!=4:
    print "python CombineRNA.py inDir outfile method"
    sys.exit()

inDir = sys.argv[1]
method = sys.argv[3]
if not os.path.exists( inDir ):
    print inDir+" not found"
    sys.exit()

if method =="RSEM_Hugo":
    RSEM_Hugo_TOIL_norm_counts (inDir, sys.argv[2])
elif method =="Kallisto_est_counts":
    Kallisto_est_counts (inDir,sys.argv[2])
elif method =="Kallisto_tpm":
    Kallisto_tpm (inDir,sys.argv[2])
elif method =="RSEM_tpm":
    RSEM_tpm (inDir,sys.argv[2])
elif method =="RSEM_IsoPct":
    RSEM_IsoPct (inDir,sys.argv[2])
else:
    print "bad method"
    sys.exit()
