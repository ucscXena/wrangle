import string, os, sys

def findTarball(tarball):
    if tarball[-7:]!=".tar.gz":
        return False
    if tarball[-13:]=="wiggle.tar.gz":
        return False
    return True

def buildLocalDir (dirname,method):
    if not os.path.exists(dirname):
        os.system("mkdir "+dirname)
    items = string.split(method,'/')
    if len(items)==1:
        if not os.path.exists(dirname+"/"+method):
            os.system("mkdir "+dirname+"/"+method)
    elif len(items)==2:
        if not os.path.exists(dirname+"/"+items[0]):
            os.system("mkdir "+dirname+"/"+items[0])
        if not os.path.exists(dirname+"/"+items[0]+"/"+items[1]):
            os.system("mkdir "+dirname+"/"+items[0]+"/"+items[1])
    else:
        print "ERROR: dir setup problem"
        sys.exit()

def getData(tarball, localDir, method, suffixDic):
    filename = string.replace(string.split(tarball,"/")[-1],".tar.gz","")
    os.system("s3cmd get --requester-pays --force "+ tarball +" "+localDir+"/")
    os.system("tar -zxf "+localDir+"/"+filename+".tar.gz -C "+localDir)
    os.system("rm "+localDir+"/"+filename+".tar.gz")
    for suffix in suffixDic[method]:
        datafile = localDir+"/"+filename+"/"+method+"/"+filename+suffix
        os.system("mv -f "+ datafile+" "+localDir+"/"+method+"/")
        print datafile
    os.system("rm -rf " +localDir+"/"+filename)

if len(sys.argv[:])!=3:
    print "python getTOILdata.py listinput localdir"
    sys.exit()

METHOD="RSEM"

SUFFIX = {
    "RSEM":[".rsem_isoforms.results",".rsem_genes.results"], #".rsem.isoform.norm_counts.tab"],
    "RSEM/Hugo":[".rsem.genes.norm_counts.hugo.tab"],
    "Kallisto":[".abundance.tsv"]
}

fin =open(sys.argv[1],'r')
localdir = sys.argv[2]
buildLocalDir(localdir, METHOD)

for line in fin.readlines():
    tarball =string.split(line)[3]
    if not findTarball (tarball):
        continue
    print tarball
    getData(tarball, localdir, METHOD, SUFFIX)
