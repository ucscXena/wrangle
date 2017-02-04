import string, sys, os
import math

def parseUnit(unit):  # working progress
    if string.find(unit,"log") != -1:
        logBase = string.strip(string.split(unit,'(')[0])[3:]
        theta = string.strip(string.split(string.split(unit,'(')[1], ')')[0])
        theta = string.strip(string.split(theta, '+')[1])
        if logBase != '' and theta != '':
            return [int(logBase), float(theta)]
    else:
        return None

def setZeroNA_log2TMPplus (value):
    #for unit = log2(tmp+0.001)
    r = parseUnit("log2(tmp+0.001)")
    logBase, theta = r

    zeroInLog = math.log(theta, logBase)

    try:
        if abs(float(value) - zeroInLog) < 0.001:
            return ' '
        else:
            return value
    except:
        return ' '


def process (infile, outfile, action):
    fout=open(outfile,'w')

    #header
    fin=open(infile,'rU')
    fout.write(fin.readline())

    while 1:
        line = fin.readline()
        if line == '':
            break
        data = string.split(line[:-1], '\t')
        values = map(lambda x: action(x), data[1:])
        fout.write(data[0] + '\t' + string.join(map(lambda x: str(x), values), '\t') +'\n')

    fout.close()

if len(sys.argv[:])!= 3:
    print "python transformGenomicMatrix fin fout"
    sys.exit()

infile = sys.argv[1]
outfile= sys.argv[2]

process(infile, outfile, setZeroNA_log2TMPplus)
