import string, sys

def findID_GenomicMatrix_lowdata  (inputFile, outputFile, lowData_threshold):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    #header
    line = fin.readline()

    #nCOL
    nCOL = len(string.split(line[:-1],'\t'))

    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line[:-1],'\t')
        if nCOL != len(data):
            print "Wrong nCOL", data[0]
            sys.exit()

        found = 0
        for i in range(1, nCOL):
            if string.strip(data[i]) in ['','NA'] :
                found = found + 1

        if 1 - (found / float(nCOL)) < lowData_threshold:
            fout.write(data[0]+'\n')

    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python findID_GenomicMatrix_lowdata.py input lowData_threshold(0.1) lowDataIDoutput"
    print
    sys.exit()

inputfile = sys.argv[1]
lowData_threshold = float(sys.argv[2])
outputfile = sys.argv[3]

findID_GenomicMatrix_lowdata (inputfile, outputfile, lowData_threshold)
