import string,os,sys

if len(sys.argv[:])!=3:
    print "python download.py listfile dir"
    sys.exit()

flist = open(sys.argv[1],'r')
targetDir  = sys.argv[2]+"/"
for line in flist.readlines():
    repo, id = string.split(line)
    if string.find(repo,"tcga")!=-1:
        TCGA =1
        key = "HAUSSL_bionimbus_pcawg_tcga.pem"
    else:
        TCGA =0
        key = "icgckeyfile.txt"
    os.system("./cghub/bin/gtdownload -c "+ key +" -d "+ repo+"cghub/data/analysis/download/"+id)
    os.system("gunzip "+id+"/*vcf.gz")
    os.system("mv "+id+"/*vcf " +targetDir)
    os.system("rm -rf "+id)
    os.system("rm "+ id+".gto")
