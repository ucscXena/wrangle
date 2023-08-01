import sys, math

def Log2xplusTheta  (inputFile, outputFile,theta):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    #header
    line = fin.readline()
    fout.write(line)

    #nCOL
    nCOL = len(line[:-1].split('\t'))

    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = line[:-1].split('\t')
        if nCOL != len(data):
            print ("Wrong nCOL", data[0])
            sys.exit()

        fout.write(data[0])

        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                fout.write("\tNA")
            else:
                value = math.log((float(data[i]) + theta),2)
                value = "%.4f" % (value)
                fout.write("\t"+value)
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print ("python Log2xplusTheta.py input(nolog) output(log2(x+theta)) theta")
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
theta = float(sys.argv[3])
Log2xplusTheta (inputfile, outputfile, theta)
