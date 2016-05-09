import string, os, sys,json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")
import CGDataLib

ASSEMBLY ="hg19"
DNA_VAF_cutoff = 0.00

def findCandidateInputDatasets (dataType, dataDirs):
    datasets=[]
    for indir in dataDirs:
        for root,dir,files in os.walk(indir,followlinks=True):
            if root[-1]=="/":
                root =root[:-1]
            if files ==[]:
                continue
            for file in files:
                if file== "cohort.json":
                    J = json.loads(open(root+"/"+file).read())
                    name = J["name"]
                    if J.has_key("default"):
                        if J["default"].has_key(dataType):
                            default =J["default"][dataType]
                            if default !=[]:
                                for ds in default:
                                    if ds not in datasets:
                                        datasets.append(ds)

    return datasets

def findQulifiedMutationVectorDatasets (bookDic, candidateDatasets):
    mutationVectorDatasets=[]
    genomicMatrixDatasets=[]

    for name in bookDic.keys():
        if name in candidateDatasets:
            #print name
            path = bookDic[name]["path"]
            type = bookDic[name]["type"]

            if not os.path.exists(path):
                continue
            if not os.path.exists(path+".json"):
                print "no .json"
                continue

            if type !="mutationVector" :
                continue

            J = json.loads(open(path+".json").read())
            if not J.has_key("assembly") or  J["assembly"]!=ASSEMBLY:
                print name, "no assembly or wrong assebly"
            else:
                mutationVectorDatasets.append(path)

            if not os.path.exists(path+"_gene") or not os.path.exists(path+"_gene.json"):
                print name, "no gene level file or no gene level json file"
            else:
                genomicMatrixDatasets.append(path+"_gene")
    
    return mutationVectorDatasets, genomicMatrixDatasets


def outputMutationVector (outfile, mutationVectorDatasets):
    HEADER=""
    DNA_VAF_col=-1
    for path in mutationVectorDatasets:
        fin = open(path,'r')
        header = fin.readline()
        if header[0]=="#":
            header = header[1:]
        if HEADER=="":
            HEADER= header
            data = string.split(header,'\t')
            for i in range(0,len(data)):
                if data[i] =="DNA_VAF":
                    DNA_VAF_col =i
            fout=open(outfile,'w')
            fout.write(HEADER)
        if header != HEADER:
            print "error, different header\n"
            continue
        while 1:
            line = fin.readline()
            if line =="":
                break
            data = string.split(line,'\t')
            DNA_VAF = data[DNA_VAF_col]
            if DNA_VAF!="":
                if float(DNA_VAF) > DNA_VAF_cutoff:
                    fout.write(line)
            else:
                fout.write(line)
        fin.close()
    fout.close()


if (len(sys.argv[:]) < 4):
    print "python virtualData.py dataType(e.g.somatic mutations) outputfile dataDir(s)"
    sys.exit()

dataType = sys.argv[1]
if dataType not in ["somatic mutation"]:
    print "unknown dataType"
    sys.exit()


outfile = sys.argv[2]
dataDirs = sys.argv[3:]

candidateDatasets = findCandidateInputDatasets (dataType, dataDirs)

bookDic={}
for dir in dataDirs:
    dic = CGDataLib.cgWalk(dir,1)
    bookDic.update(dic)

#print dataType
#print candidateDatasets

mutationVectorDatasets, genomicMatrixDatasets = findQulifiedMutationVectorDatasets (bookDic, candidateDatasets)

print len(mutationVectorDatasets)
print len(genomicMatrixDatasets)

if len(mutationVectorDatasets)!=0:
    outputMutationVector (outfile, mutationVectorDatasets)
    
if len(genomicMatrixDatasets)!=0:
    os.system("python mergeGenomicMatrixFiles_memEfficient.py "+outfile+"_gene . " +string.join(genomicMatrixDatasets," "))

