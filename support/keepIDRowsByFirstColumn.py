import string, sys,os

def listing (mfile):
    fin = open(mfile,'r')
    dic= {}
    for line in fin.readlines():
        data = string.strip(string.split(line,'\t')[0])
        if data not in dic:
            dic[data]=0
    fin.close()
    return dic

def keepProbes (inputFile, outputFile, keep_dic):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    line = fin.readline()
    fout.write(line)
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line,'\t')
        if data[0] not in keep_dic:
            continue
        else:
            fout.write(line)
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python keepIDRowsByFirstColumn.py infile output remove_list(first_column_id)"
    sys.exit()

listfile = sys.argv[3]
keep_dic  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

keepProbes (inputfile, outputfile, keep_dic)
