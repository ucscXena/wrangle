import string, sys, os

def cleanProbeMapByData (matrixfile, probeMapfile, outputfile):
    #matrix features
    features = {}
    fin = open(matrixfile,'r')
    while 1:
        line = fin.readline()
        if line == "":
            break
        data = string.split(line[:-1],'\t')
        f = data[0]
        #print f
        for item in data[1:]:
            if item not in ["","NA","0"]:
                #print item
                features[f]=''
                break
    fin.close()

    fprobe = open(probeMapfile,'U')
    fout = open(outputfile,'w')

    #probeMap
    line = fprobe.readline()
    fout.write(line)

    while 1:
        line = fprobe.readline()
        if line =="":
            break

        data = string.split(line[:-1],'\t')
        if data[0] in features:
           fout.write(line)

    fprobe.close()
    fout.close()

if len(sys.argv[:])!=4:
    print "python cleanProbeMapByData.py genomicMatrix probeMap newProbeMap"
    print
    sys.exit()

matrixfile = sys.argv[1]
probeMapfile = sys.argv[2]
outputfile = sys.argv[3]

cleanProbeMapByData (matrixfile, probeMapfile, outputfile)
