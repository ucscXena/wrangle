import string, sys, os

def cleanProbeMapByData (matrixfile, probeMapfile, outputfile):
    fprobe = open(probeMapfile,'U')
    fout = open(outputfile,'w')

    #matrix features
    os.system("cut -f 1 "+ matrixfile + " > features")
    fin = open("features",'r')
    features = {}
    for f in fin.readlines():
        f = string.strip(f)
        features[f]=''
    fin.close()

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
