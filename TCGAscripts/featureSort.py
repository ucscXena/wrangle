import string, pprint,sys
from TCGAfeature import *

def longTitle(new,change):
    index= featureLongTitle.keys()
    index.sort()
    for key in new.keys():
        if key not in index:
            if change.has_key(key):
                featureLongTitle[key]=change[key]
            else:
                featureLongTitle[key]=new[key]

            
def shortTitle(new,change):
    index= featureShortTitle.keys()
    index.sort()
    for key in new.keys():
        if key not in index:
            if change.has_key(key):
                featureShortTitle[key]=change[key]
            else:
                featureShortTitle[key]=new[key]

fin=open("missing",'r')
new={}
for line in fin.readlines():
    line=string.strip(line)
    data =string.split(line,"\t")
    new[data[1]]=data[2]
fin.close()

change={}
fin = open("change",'r')
for line in fin.readlines():
    line=string.strip(line)[34:]
    data =string.split(line," new name: ")
    change[data[1]]=data[0]
fin.close()

longTitle(new,change)
shortTitle(new,change)

fout=open("missingpython",'w')
sys.stdout = fout
fout.write("featureLongTitle=\n")
pprint.pprint(featureLongTitle)

fout.write("\n")
fout.write("featureShortTitle=\n")
pprint.pprint(featureShortTitle)

fout.close()

sys.stdout = sys.__stdout__ 
