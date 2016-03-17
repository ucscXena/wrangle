import string,os,sys

freezefile = "release_mar2016.v1.txt"
keyword ="_variant_calling_repo"

def parseFreezeFile (infile):
    fin =open (infile,'r')
    reg ={}
    
    #header
    data = string.split(fin.readline(),'\t')
    N = len(data)
    for i in range (0, N):
        if string.find(data[i],keyword)!=-1:
            method = string.replace(data[i],keyword,"")
            if not os.path.exists("PCAWG/"+method):
                os.system("mkdir PCAWG/"+ method)
            reg[i]=open(method+"_list",'w')
            

    #datalines
    for line in fin.readlines():
        data = string.split(line,'\t')
        n =len(data)
        for i in range (0, N):
            if i in reg:
                if string.strip(data[i])!="":
                    fout = reg[i]
                    fout.write(data[i]+'\t'+data[i+1]+'\t'+data[i+2]+'\n')

    fin.close()
    for i in reg:
        fout= reg[i]
        fout.close()

parseFreezeFile(freezefile)
