import string
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomicMatrixInput", type=str)
    parser.add_argument("probeMap", type=str)
    parser.add_argument("geneMatrixOutput", type=str)
    args = parser.parse_args()

    probeMap = readProbeMap(args.probeMap)
    geneMatrix, colN = process (args.genomicMatrixInput, probeMap, args.geneMatrixOutput)
    outputMatrix (geneMatrix, args.geneMatrixOutput, args.type, args.log2, args.theta, colN)

def process(genomicMatrixInput, probeMap, output):
    genomicFp = open(genomicMatrixInput)
    fout = open(output, 'w')

    line = genomicFp.readline()
    fout.write(line)
    
    while 1:
        line = genomicFp.readline()
        if line =='':
            break
        
        firstEnd = line.find('\t')
        probe = data[0:firstEnd]
        if probe not in probeMap:
            gene = probe
            print probe, "not found in probeMap"
        else:
            gene = probeMap[probe]
    
        fout.write(gene+ line[firstEnd:])
    genomicFp.close()
    fout.close()


def readProbeMap(probeMapFile):
    probeMap = dict()
    fp = open(probeMapFile)
    fp.readline()
    for row in fp:
        tokens = row.rstrip().split("\t")
        probeMap[tokens[0]] = string.split(tokens[1],',')[0]
    fp.close()
    return(probeMap)

if __name__ == "__main__":
    main()
