import string, sys, math

def Log2xplusThetaToBeforeLog  (inputFile, outputFile, theta):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    #header
    line = fin.readline()
    fout.write(line)

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

        fout.write(data[0])

        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                fout.write("\tNA")
            else:
                value = math.pow(2,float(data[i])) - theta
                # RNAseq specific: before log, the lowest value should be zero
                if value < theta:
                    value = 0
                value = "%.4f" % (value)
                fout.write("\t"+value)
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:]) != 4:
    print "Log2xplusThetaToBeforeLog.py input(log2(x+theta) output(beforeLog) theta"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
theta = float(sys.argv[3])

Log2xplusThetaToBeforeLog (inputfile, outputfile, theta)
