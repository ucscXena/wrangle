import os, sys, string

def process (inputdir):
    for file in os.listdir(inputdir):
        if (not os.path.isfile(file)):
            continue
        input = inputdir + "/" + file
        sample = string.replace(file,'.txt','')

        fout = open("new",'w')
        fout.write("sample\t"+ sample+'\n')
        fout.close()
        os.system("cat "+ input + " >> new") 
        print sample
        os.system("mv new "+ inputdir + "/" + sample +".xena")

if len(sys.argv[:]) != 2:
    print "python DNAmethyl.py inputdir"
    sys.exit()

inputdir = sys.argv[1]
if inputdir[-1]=='/':
    inputdir = inputdir[:-1]
process (inputdir)
