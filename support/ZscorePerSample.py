import string, sys, math

def ZscorePerSample  (inputFile, outputFile):
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

        # average and count
        count = 0
        total = 0.0
        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                pass
            else:
                count = count + 1
                total = total + float(data[i])
        average = total / count

        # std
        total = 0.0
        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                fout.write("\tNA")
            else:
                value = float(data[i])
                total = total + (value - average) * (value - average)
        std = math.sqrt(total /(count -1))

        fout.write(data[0])
        for i in range(1, nCOL):
            if data[i] in ["NA",""]:
                fout.write("\tNA")
            else:
                z = (float(data[i]) - average) / std
                z = "%.4f" % (value)
                fout.write("\t"+value)
        fout.write("\n")

    fin.close()
    fout.close()

if len(sys.argv[:])!=3:
    print "python ZscorePerSample.py input(genomicMatrix_nomralDistributionBest) output(zscore)"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
ZscorePerSample (inputfile, outputfile)