import sys, string
import numpy

def averageDup (inputfile, outputfile):
    fin = open(inputfile, 'U')
    fout = open(outputfile, 'w')

    #header - find dups
    dupHeader = {}
    DUP = 0
    line = fin.readline()
    data = string.split(line[:-1], '\t')
    for i in range (1, len(data)):
        sample = data[i]
        if sample not in dupHeader:
            dupHeader[sample]=[i]
        else:
            DUP =1
            dupHeader[sample].append(i)
            print sample, dupHeader[sample]

    #no dup data
    if not DUP:
        fout.write(line)
        fout.write(fin.read())
    else:
        #dup data
        samples = dupHeader.keys()
        fout.write('sample\t' + string.join(samples,'\t')+'\n')

        while 1:
            line = fin.readline()
            if line == '':
                break
            data = string.split(line[:-1], '\t')
            gene = data[0]
            fout.write(gene)

            for sample in samples:
                pos_list = dupHeader[sample]

                if len(pos_list) >1:
                    values = map(lambda p: data[p],  pos_list)
                    values = map(float, filter(lambda v: v!='' and v!= 'NA', values))
                    if len(values) == 0:
                        fout.write('\t')
                    else:
                        ave = numpy.average(values)
                        fout.write('\t' + str(ave))
                else:
                    fout.write('\t' + data[pos_list[0]])
            fout.write("\n")

    fin.close()
    fout.close()

if len(sys.argv[:]) != 3:
    print "python duplicateAverageGenomicMatrix.py matrixinput outputfile"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]

averageDup (inputfile, outputfile)
