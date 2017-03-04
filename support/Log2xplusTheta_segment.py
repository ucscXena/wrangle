import string, sys, math

def Log2xplusTheta  (inputFile, outputFile,theta):
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

        fout.write(string.join(data[0:5],'\t'))

        if data[5] in ["NA",""]:
            fout.write("\tNA")
        else:
            value = math.log((float(data[5]) + theta),2)
            value = "%.4f" % (value)
            fout.write("\t"+value)
        fout.write("\n")
    fin.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python Log2xplusTheta_segment.py input(nolog) output(log2(x+theta)) theta"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
theta = float(sys.argv[3])
Log2xplusTheta (inputfile, outputfile, theta)
