import string, os, sys,json

os.sys.path.insert(0, os.path.dirname(__file__)+"../CGDataNew")
import CGDataLib

ASSEMBLY ="hg19"

if (len(sys.argv[:]) < 4):
    print "python virtualData.py dataType outputfile dir(s)"
    sys.exit()

dataType = sys.argv[1]
if dataType not in ["somatic mutation"]:
    print "unknown dataType"
    sys.exit()

print dataType
bookDic={}
for dir in sys.argv[3:]:
    dic = CGDataLib.cgWalk(dir,1)
    bookDic.update(dic)

datasets=[]
for indir in sys.argv[3:]:
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
HEADER=""
for name in bookDic.keys():
    if name in datasets:
        print name
        path = bookDic[name]["path"]
        type = bookDic[name]["type"]
        if not os.path.exists(path):
            continue
        if dataType in ["somatic mutation"]:
            if type !="mutationVector":
                continue
        if not os.path.exists(path+".json"):
            print "no .json"
            continue
        J = json.loads(open(path+".json").read())
        if not J.has_key("assembly"):
            print "no assembly"
            continue
        if J["assembly"]!=ASSEMBLY:
            print "wrong assembly"
            continue

        if dataType in  ["somatic mutation"]:
            fin = open(path,'r')
            header = fin.readline()
            if header[0]=="#":
                header = header[1:]
            if HEADER=="":
                HEADER= header
                fout=open(sys.argv[2],'w')
                fout.write(HEADER)
                fout.close()
            if header != HEADER:
                print "error, different header\n"
                continue
            if dataType =="somatic mutation":
                os.system ("more +2 "+ path +" >> "+ sys.argv[2])
