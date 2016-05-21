import string,sys, json, os

if len(sys.argv[:])!=2:
    print "python jsonLint.py dirToLint"
    sys.exit()

inDir = sys.argv[1]
print inDir
for root,dir,files in os.walk(inDir,followlinks=True):
    if root[-1]=="/":
        root =root[:-1]
    if files ==[]:
        continue
    for file in files:
        if file[-5:]!=".json":
            continue
        print file
        J = json.loads(open(root+"/"+file).read()) 
