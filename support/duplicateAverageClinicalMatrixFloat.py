import sys, string, os
import numpy

def detectDup (inputfile):
    r = os.popen("cut -f 1 " + inputfile + " | sort | uniq -c | sort -k 1 | grep -v \' \'1\' \' |wc -l ", 'r')
    dupNum = int(string.strip(r.read()))
    return dupNum

def averageDup (inputfile, outputfile):
    fin = open(inputfile, 'U')
    fout = open(outputfile, 'w')

    #header
    line = fin.readline()
    fout.write(line)
    headers = string.split(line[:-1],'\t')

    dic = {}
    samples =[]
    #read the whole file into memory
    for line in fin.readlines():
        data = string.split(line, '\t')
        sample = data[0]
        if sample not in dic:
            dic[sample] ={}
            samples.append(sample)
        for i in range(1, len(headers)):
            header = headers[i]
            if header not in dic[sample]:
                dic[sample][header] = []
            data[i] = string.strip(data[i])
            if data[i] not in dic[sample][header]:
                dic[sample][header].append(data[i])

    #output resolve dup
    for sample in samples:
        fout.write(sample)
        for header in headers[1:]:
            values = dic[sample][header]
            for v in values:
                if v in ['', 'NA']:
                    values.remove(v)
            if len(values) == 0:
                value = ''
            elif len(values) == 1:
                value = values[0]
            else:
                allFloats = 1
                for v in values:
                    try:
                        v = float(v)
                    except:
                        value = string.join(values,',')
                        allFloats = 0
                if allFloats:
                    values = map(float, values)
                    ave = numpy.average(values)
                    value = str(ave)
            fout.write('\t' + value)
        fout.write("\n")

    fin.close()
    fout.close()

if len(sys.argv[:]) != 3:
    print "python duplicateAverageClinicalMatrixFloat.py clinMatrixinput outputfile"
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]

if detectDup (inputfile):
    print "there is dups"
else:
    print "no dups"
    sys.exit()

averageDup (inputfile, outputfile)
