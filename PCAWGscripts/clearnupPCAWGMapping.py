import string,os,sys

if len(sys.argv[:])!=3:
    print "python clearnupPCAWGMapping.py input outputClean"
    sys.exit()

fin = open(sys.argv[1],'r')
fout = open(sys.argv[2],'w')

fout.write(fin.readline())
for line in fin.readlines():
    print line
    olds, news = string.split(string.strip(line),'\t')
    olds = string.replace(olds,'"','')
    news = string.replace(news,'"','')
    olds = string.split(olds,',')
    news = string.split(news,',')
    if (len(olds) != len(news)):
        if len(olds)==1:
            print olds, news
            for i in range(0,len(news)):
                fout.write(olds[0]+'\t'+news[i]+'\n')
        else:    
            print olds, news
        continue
    else:
        for i in range(0,len(olds)):
            fout.write(olds[i]+'\t'+news[i]+'\n')
fin.close()
fout.close()
