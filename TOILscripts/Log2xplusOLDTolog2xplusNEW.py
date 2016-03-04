import string, sys, math

def Log2xplusOLDTolog2xplusNEW  (inputFile, outputFile, oldtheta, newtheta):
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
                x = math.pow(2,float(data[i]))- oldtheta
                value = math.log((x + newtheta),2)
                value = "%.4f" % (value)
                fout.write("\t"+value)
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:])!=5:
    print "Log2xplusOLDTolog2xplusNEW.py input(log2(x+old) output(log2(x+new)) oldtheta newtheta"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
oldtheta = float(sys.argv[3])
newtheta = float(sys.argv[4])
Log2xplusOLDTolog2xplusNEW (inputfile, outputfile, oldtheta, newtheta)
