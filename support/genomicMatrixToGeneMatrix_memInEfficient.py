import string
import argparse
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomicMatrixInput", type=str)
    parser.add_argument("probeMap", type=str)
    parser.add_argument("geneMatrixOutput", type=str)
    parser.add_argument("type", type=str, choices=('add', 'average'), help='add or average')
    parser.add_argument("log2", type=str, choices=('log2', 'nolog'), help='log or nolog')
    parser.add_argument("theta", type=float, default=0, help='theta value')
    args = parser.parse_args()

    probeMap = readProbeMap(args.probeMap)
    geneMatrix, colN = process (args.genomicMatrixInput, probeMap, args.type, args.geneMatrixOutput)
    outputMatrix (geneMatrix, args.geneMatrixOutput, args.type, args.log2, args.theta, colN)

def process(genomicMatrixInput, probeMap, type, output):
    geneMatrix = {}
    genomicFp = open(genomicMatrixInput)

    fout = open(output, 'w')
    line = genomicFp.readline()
    N = len(string.split(line,'\t'))
    fout.write(line)
    fout.close()

    while 1:
        line = genomicFp.readline()
        if line =='':
            break
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


def ADD_type_conversion (list, log2, theta, colN):
    ret_list = []
    for i in range(0, colN):
        ret_list.append('NA')

    for row in list:
        for i in range(0, colN):
            try:
                value =  float(row[i])
                if log2:
                    value = math.pow(2, value) - theta
            except:
                continue
            if ret_list[i] == 'NA':
                ret_list[i] = value
            ret_list[i] = ret_list[i] + value
    if log2:
        for i in range(0, colN):
            if ret_list[i]!= 'NA':
                ret_list[i] = math.log(ret_list[i] + theta, 2)
    return ret_list

def Average_type_conversion (list, log2, theta, colN):
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
                if log2:
                    value = math.pow(2, value) - theta
            except:
                continue
            if total_list[i] == 'NA':
                total_list[i] = value
            total_list[i] = total_list[i] + value
            count_list[i] = count_list[i] + 1

    for i in range(0, colN):
        if total_list[i]!='NA':
            ret_list [i] = total_list[i] / count_list [i]

    if log2:
        for i in range(0, colN):
            if ret_list[i]!= 'NA':
                ret_list[i] = math.log(ret_list[i] + theta, 2)

    return ret_list


def outputMatrix (geneMatrix, output, type, log2, theta, colN):
    fout = open(output, 'a')
    for gene in geneMatrix:
        geneData = geneMatrix[gene]
        if type == 'add':
            ret_list = ADD_type_conversion (geneData, log2, theta, colN)
        elif type == 'average':
            ret_list = Average_type_conversion (geneData, log2, theta, colN)
        fout.write(gene)
        for i in range(0, colN):
            fout.write('\t'+str(ret_list[i]))
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
