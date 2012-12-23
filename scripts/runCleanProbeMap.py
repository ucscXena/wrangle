import string, os, sys, json
import optparse

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")

parser = optparse.OptionParser()
parser.add_option("--inDir", action="store", type="string", dest="inDir")
parser.add_option("--outDir", action="store", type="string", dest="outDir")
parser.add_option("--file", action="store", type="string", dest="infile")
(options, args) = parser.parse_args()


#print options

def printUsage():
    print "python runFlatten.py --inDir=inputDir --outDir=outputDir\n"
    print "options:"
    print "          --file=input probemap file"


def process(dirpath, file, outDir):
    fin = open(dirpath+file,'r')
    lines = fin.readlines()
    WRITE=0
    for i in range(len(lines)-1, -1, -1):
        EDIT=0
        chrom = string.split(lines[i],'\t')[2]
        if chrom in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                     "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                     "chr22","chrX","chrY"]:
            continue
        if chrom[0:3]=="CHR":
            chrom="chr"+chrom[3:]
            EDIT =1
        if chrom =="chr23":
            chrom= "chrX"
            EDIT =1
        if chrom =="chr24":
            chrom= "chrY"
            EDIT =1
        if chrom in ["1","01"]:
            chrom ="chr1"
            EDIT =1
        if chrom in ["2","02"]:
            chrom ="chr2"
            EDIT =1
        if chrom in ["3","03"]:
            chrom ="chr3"
            EDIT =1
        if chrom in ["4","04"]:
            chrom ="chr4"
            EDIT =1
        if chrom in ["5","05"]:
            chrom ="chr5"
            EDIT =1
        if chrom in ["6","06"]:
            chrom ="chr6"
            EDIT =1
        if chrom in ["7","07"]:
            chrom ="chr7"
            EDIT =1
        if chrom in ["8","08"]:
            chrom ="chr8"
            EDIT =1
        if chrom in ["9","09"]:
            chrom ="chr9"
            EDIT =1
        if chrom =="10":
            chrom ="chr10"
            EDIT =1
        if chrom =="11":
            chrom ="chr11"
            EDIT =1
        if chrom =="12":
            chrom ="chr12"
            EDIT =1
        if chrom =="13":
            chrom ="chr13"
            EDIT =1
        if chrom =="14":
            chrom ="chr14"
            EDIT =1
        if chrom =="15":
            chrom ="chr15"
            EDIT =1
        if chrom =="16":
            chrom ="chr16"
            EDIT =1
        if chrom =="17":
            chrom ="chr17"
            EDIT =1
        if chrom =="18":
            chrom ="chr18"
            EDIT =1
        if chrom =="19":
            chrom ="chr19"
            EDIT =1
        if chrom =="20":
            chrom ="chr20"
            EDIT =1
        if chrom =="21":
            chrom ="chr21"
            EDIT =1
        if chrom =="22":
            chrom ="chr22"
            EDIT =1
        if chrom =="23":
            chrom ="chrX"
            EDIT =1
        if chrom =="24":
            chrom ="chrY"
            EDIT =1
        if EDIT:
            data = string.split(lines[i],'\t')
            data[2]=chrom
            lines[i]=string.join(data,'\t')
        else:
            lines.pop(i)
        WRITE=1
    fin.close()
    if not WRITE:
        os.system("cp "+dirpath+file+" "+outDir)
    else:
        fout =open(outDir+file,"w")
        fout.writelines(lines)
        fout.close()
    os.system("cp "+dirpath+file+".json "+outDir)
    
        

if options.inDir==None or options.outDir ==None:
    printUsage()
    sys.exit()

inDir = options.inDir
outDir =options.outDir

if inDir[-1]!="/":
    inDir = inDir +"/"

if outDir[-1]!="/":
    outDir = outDir +"/"
    
    
if inDir == outDir:
    print "inDir and outDir is the same"
    sys.exit()
    
#inDir="data/probeMap/"
#outDir="data_flatten/probeMap/"


for dirpath, dirnames, filenames in os.walk(inDir, followlinks=True):
    for file in filenames:
        if os.path.abspath(dirpath)==os.path.abspath(outDir):
            continue
        if options.infile != None and options.infile != file:
            continue
        path=dirpath+file
        if not os.path.exists(path+".json"):
            print path,"missing .json"
            continue

        J = json.loads(open(path+".json").read())
        if J["type"]=="probeMap":
            print path
            process(dirpath, file, outDir)

