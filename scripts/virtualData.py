import string, os, sys,json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")
import CGDataLib

ASSEMBLY ="hg19"
DNA_VAF_cutoff = 0.00

if (len(sys.argv[:]) < 5):
    print "python virtualData.py dataType(e.g.somatic mutations) outFileType(e.g.mutationVector) outputfile dir(s)"
    sys.exit()

dataType = sys.argv[1]
if dataType not in ["somatic mutation"]:
    print "unknown dataType"
    sys.exit()

outputFileType = sys.argv[2]

print dataType
bookDic={}
for dir in sys.argv[4:]:
    dic = CGDataLib.cgWalk(dir,1)
    bookDic.update(dic)

datasets=[]
for indir in sys.argv[4:]:
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

print datasets
genomicMatrixDatasets=[]
mutationVectorDatasets=[]

for name in bookDic.keys():
    if name in datasets:
        print name
        path = bookDic[name]["path"]
        type = bookDic[name]["type"]

        if not os.path.exists(path):
            continue
        if not os.path.exists(path+".json"):
            print "no .json"
            continue

        if type == outputFileType :
            pass
        elif type =="mutationVector" and outputFileType =="genomicMatrix":
            pass
        else:
            continue

        J = json.loads(open(path+".json").read())

        if dataType in  ["somatic mutation"] and outputFileType == "mutationVector":
            if not J.has_key("assembly"):
                print "no assembly"
                continue
            if J["assembly"]!=ASSEMBLY:
                print "wrong assembly"
                continue
            mutationVectorDatasets.append(path)

        if dataType in  ["somatic mutation"] and outputFileType == "genomicMatrix":
            genomicMatrixDatasets.append(path+"_gene")

if len(mutationVectorDatasets)!=0:
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
            fout=open(sys.argv[3],'w')
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
            #os.system ("more +2 "+ path +" >> "+ sys.argv[3])
        fin.close()
    fout.close()

if len(genomicMatrixDatasets)!=0:
    os.system("python mergeGenomicMatrixFiles.py "+ sys.argv[3]+" " + string.join(genomicMatrixDatasets," "))
