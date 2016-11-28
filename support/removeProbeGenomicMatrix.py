import string, sys,os

def listing (mfile):
    fin = open(mfile,'r')
    list =[]
    for line in fin.readlines():
        data = string.strip(line)
        if data not in list:
            list.append(data)
    fin.close()
    return list

def removeProbes (inputFile, outputFile, remove_list):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')
    line = fin.readline()
    fout.write(line)
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line,'\t')
        if data[0] in remove_list:
            continue
        else:
            fout.write(line)
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python removeRowsByFirstColumn.py infile output remove_list(one_id_per_line)"
    sys.exit()

listfile = sys.argv[3]
remove_list  = listing (listfile)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

removeProbes (inputfile, outputfile, remove_list)
