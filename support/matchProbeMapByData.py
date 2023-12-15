import sys, os

def matchProbeMapByData (matrixfile, probeMapfile):
    #matrix features
    features = {}
    fin = open(matrixfile,'r')
    while 1:
        line = fin.readline()
        if line == "":
            break
        data = line[:-1].split('\t')
        feature = data[0]
        features[f]=''
    fin.close()
    matrix_set = {features.keys()}
    
    #probeMap
    features = {}
    fprobe = open(probeMapfile,'U')
    line = fprobe.readline()
    while 1:
        line = fprobe.readline()
        if line =="":
            break

        data = line[:-1].split('\t')
        feature = data[0]
        features[f]=''
    fprobe.close()
    probemap_set = {features.keys()}

    # overlap
    overlap_set = set.intersection(matrix_set, probemap_set)

    # stats
    print ("overlap:", len(overlap_set))
    print ("matrix:", len(matrix_set), len(overlap_set)/len(matrix_set) )
    print ("probemap:", len(probemap_set), len(probemap_set)/len(matrix_set))


if len(sys.argv[:])!=3:
    print 
    print ("python matchProbeMapByData.py genomicMatrix probeMap\n")
    sys.exit()

matrixfile = sys.argv[1]
probeMapfile = sys.argv[2]

mathProbeMapByData (matrixfile, probeMapfile)
