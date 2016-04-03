import string, sys,os

def copyOriginalFile (filename, outputfile):
    os.system("cp " +filename +" "+ outputfile)

def listing (mfile):
    fin = open(mfile,'r')
    list =[]
    for line in fin.readlines():
        data = string.strip(line)
        if data not in list:
            list.append(data)
    fin.close()
    return list

def removeIDs (inputFile, outputFile, remove_list):
    fin = open(inputFile,'U')
    
    badPos=[]
    ids = string.split(fin.readline()[:-1],'\t')
    for i in range (0, len(ids)):
        for badid in remove_list:
            if string.find(ids[i], badid)!=-1:
                badPos.append(i)
                print ids[i],badid, i
                break
    fin.close()

    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line[:-1],'\t')
        fout.write(data[0])
        for i in range (1, len(data)):
            if i not in badPos:
                fout.write("\t"+data[i])
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "xenaIDreplace.py matrixfile output remove_list(one_id_per_line)"
    sys.exit()

listfile = sys.argv[3]
remove_list  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

removeIDs (inputfile, outputfile, remove_list)
