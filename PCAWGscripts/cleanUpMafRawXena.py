import string, sys
import numpy

#cleaned up pos=7 (zero) 8th column VAF

def process(inputfile, outputfile):
    fin = open(inputfile,'r')
    fout = open(outputfile,'w')
    fout.write(fin.readline())

    while 1:
        line = fin.readline()
        if line == "":
            break
        data = string.split(line,'\t')
        if string.find(data[7],'|') != -1:
            average = numpy.average(map(float, string.split(data[7])))
            data[7]= str(average)
            fout.write(string.join(data,'\t'))
        else:
            fout.write(line)

    fin.close()
    fout.close()

if __name__ == '__main__' and len(sys.argv[:])!=3:
    print "python cleanUpMafRawXena.py xenaSVMutationVectorFile(derived from maf) cleanedup_outputfile"
    sys.exit()

input = sys.argv[1]
output = sys.argv[2]

process(input, output)
