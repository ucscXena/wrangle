import sys, string
import math, numpy

def averageDup (inputfile, outputfile, LOG, theta):
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
        os.system("cp " + inputfile + " " + outputfile)
    #dup data
    else:
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
                        if LOG:
                            values = map( lambda x: math.pow(x,2) - theta, values)
                        ave = numpy.average(values)
                        if LOG:
                            ave = math.log( (ave + theta), 2)
                        fout.write('\t' + str(ave))
                else:
                    fout.write('\t' + data[pos_list[0]])
            fout.write("\n")

    fin.close()
    fout.close()

if len(sys.argv[:]) not in [3, 5]:
    print "python duplicateAverageGenomicMatrix.py matrixinput outputfile LOG2(0,1) Theta"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
LOG = 0
theta =0

if len(sys.argv[:]) == 5:
    LOG = sys.argv[3]
    if LOG not in ["0", "1"]:
        print "LOG must be 0 or 1\n"
        sys.exit()
    LOG = int(LOG)
    theta = float(sys.argv[4])

averageDup (inputfile, outputfile, LOG, theta)
