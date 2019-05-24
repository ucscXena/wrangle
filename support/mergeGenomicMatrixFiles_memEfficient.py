import sys,string,os
import json,datetime
import math
import copy
import uuid

def merge (infile_list, outfile, root ="./"):
    genes={} #key: gene value position 0,1,2, ...
    allgenes={} #key: position 0,1,2, ... value: gene
    files=[]

    if root[-1]!="/":
        root = root +"/"
    tmpDir =root + str(uuid.uuid4())

    os.system("mkdir "+ tmpDir)

    tmpList =[]
    for infile in infile_list:
        if not os.path.exists(infile):
            print infile, "not exists, skip"
        else:
            tmpList.append(infile)

    infile_list = tmpList

    for infile in infile_list:
        genes = setGeneOrder (infile, genes, tmpDir)

    allgenes = posToGene(genes)
    tmpfile = tmpDir+"/0_id"
    outputallgenes(allgenes,tmpfile)
    files.append(tmpfile)

    c=0
    for infile in infile_list:
        c=c+1
        nCOLs = getColumnSize (infile)
        cur_genes={}
        cur_genes = MapGeneOrder (infile,cur_genes, tmpDir)

        for i in range (2,nCOLs+1, 250):
            tmpfile = tmpDir+"/"+str(c)+"_tmp_"+str(int(i/250.0))
            os.system("cut -f "+ str(i)+"-"+str(i+249)+" "+infile +" > "+tmpDir+"/tmp")
            process(cur_genes, tmpDir+"/tmp", allgenes, tmpfile)
            files.append(tmpfile)

    #paste all together
    k = 100
    os.system("paste "+ string.join(files[:k], " ") + " > " + outfile)
    for i in range(k, len(files), k):
        tmpfile = str(uuid.uuid4())
        os.system("paste "+ outfile + ' ' + string.join(files[i: i+k], " ") + " > " + tmpfile)
        os.system("mv " + tmpfile + ' ' + outfile)

    os.system("rm -rf "+ tmpDir)
    return

def posToGene(genes):
    allGenes ={}
    for gene in genes:
        p = genes[gene]
        allGenes[p]=gene
    return allGenes

def outputallgenes(allgenes,outfile):
    fout=open(outfile,'w')
    fout.write("xena_sample\n")
    for i in range(0,len(allgenes)):
        gene = allgenes[i]
        fout.write(gene+"\n")
    fout.close()

def setGeneOrder (infile,genes, tmpDir):
    os.system("cut -f 1 "+infile +" > "+ tmpDir+"/tmpid")
    fin=open(tmpDir+"/tmpid",'U')
    fin.readline()
    for line in fin.readlines():
        hugo = line[:-1]
        if hugo not in genes:
            p=len(genes)
            genes[hugo]=p
    fin.close()
    return genes

def MapGeneOrder (infile, genes, tmpDir):
    os.system("cut -f 1 "+infile +" > "+ tmpDir+"/tmpid")
    fin=open(tmpDir+"/tmpid",'U')
    fin.readline()
    lines = fin.readlines()
    for i in range (0, len(lines)):
        hugo = lines[i][:-1]
        genes[hugo]= i
    fin.close()
    return genes

def getColumnSize (infile):
    fin=open(infile,'U')
    line =fin.readline()
    fin.close()
    return len(string.split(line[:-1],"\t"))

def process(cur_genes, infile, allgenes, outfile):
    fin=open(infile,'r')
    fout=open(outfile,'w')

    line = fin.readline()
    fout.write(line)

    nCOL = len(string.split(line,'\t'))
    emptyList =[]
    for i in range(0,nCOL):
        emptyList.append("")
    emptyLine =string.join(emptyList,'\t')+"\n"

    lines = fin.readlines()

    for i in range(0,len(allgenes)):
        gene = allgenes[i]
        if gene in cur_genes:
            p = cur_genes[gene]
            fout.write(lines[p])
        else:
            fout.write(emptyLine)
    fin.close()
    fout.close()
    return

if __name__ == "__main__" and len(sys.argv[:])<4:
    print "python mergeGenomicMatrixFiles_memEfficient.py outfile tmpDirRoot infiles"
    sys.exit()

if __name__ == "__main__":
    outfile = sys.argv[1]
    root = sys.argv[2]
    infile_list =sys.argv[3:]
    print root, infile_list, outfile
    merge (infile_list, outfile, root)
