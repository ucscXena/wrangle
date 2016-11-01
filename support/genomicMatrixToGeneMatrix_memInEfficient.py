import string
import argparse
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomicMatrixInput", type=str)
    parser.add_argument("probeMap", type=str)
    parser.add_argument("geneMatrixOutput", type=str)
    parser.add_argument("type", type=str, choices=('add', 'average'), help='add or average')
    args = parser.parse_args()

    probeMap = readProbeMap(args.probeMap)
    geneMatrix, colN = process (args.genomicMatrixInput, probeMap, args.type, args.geneMatrixOutput)
    outputMatrix (geneMatrix, args.geneMatrixOutput, args.type, colN)

def process(genomicMatrixInput, probeMap, type, output):
    geneMatrix = {}
    genomicFp = open(genomicMatrixInput)

    fout = open(output, 'w')
    line = genomicFp.readline()
    N = len(string.split(line,'\t'))
    fout.write(genomicFp.readline())
    fout.close()

    while 1:
        line = genomicFp.readline()
        data = string.split(line,'\t')
        probe = data[0]
        if probe not in probeMap:
            continue
        gene = probeMap[probe]
        if gene not in geneMatrix:
            geneMatrix[gene]=[]
        geneMatrix[gene].append(data[1:])
    genomicFp.close()
    return geneMatrix, N


def ADD_type_conversion (list, colN):
    ret_list = []
    for i in range(0, colN):
        ret_list.append('NA')

    for row in list:
        for i in range(0, colN):
            try:
                value =  float(row[i])
            except:
                continue
            if ret_list[i] == 'NA':
                ret_list[i] = value
            ret_list[i] = ret_list[i] + value
    return ret_list

def Average_type_conversion (list, colN):
    total_list = []
    count_list = []
    ret_list =[]

    for i in range(0, colN):
        total_list.append('NA')
        ret_list.append('NA')
        count_list.append(0)

    for row in list:
        for i in range(0, colN):
            try:
                value =  float(row[i])
            except:
                continue
            if total_list[i] == 'NA':
                total_list[i] = value
            total_list[i] = total_list[i] + value
            count_list[i] = count_list[i] + 1

    for i in range(0, colN):
        if total_list[i]!='NA':
            ret_list [i] = total_list[i] / count_list [i]

    return ret_list


def outputMatrix (geneMatrix, output, type, colN):
    fout = open(output, 'a')
    for gene in matrix:
        geneData = matrix[gene]
        if type == 'add':
            ret_list = ADD_type_conversion (list, colN)
        elif type == 'average':
            ret_list = Average_type_conversion (list, colN)
        fout.write(gene)
        for i in range(0, colN):
            fout.write('\t'+ret_list[i])
        fout.write('\n')
    fout.close()

def readProbeMap(probeMapFile):
    probeMap = dict()
    fp = open(probeMapFile)
    fp.readline()
    for row in fp:
        tokens = row.rstrip().split("\t")
        probeMap[tokens[0]] = tokens[1]
    fp.close()
    return(probeMap)

if __name__ == "__main__":
    main()
