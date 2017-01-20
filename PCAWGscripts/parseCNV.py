import string, sys, os
dataCol = 3

def writeHeader(fout, line):
    fout.write("sampleID\tchr\tstart\tend\tvalue" +'\t')
    data = string.split(line,'\t')
    fout.write(string.join(data[3:],'\t'))

def process(inDir, output):
    fout = open(output,'w')
    count = 0
    for file in os.listdir(inDir):
        fname = inDir + file
        sampleName = string.split(file,'.')[0]

        fin = open(fname,'r')
        if count == 0:
            writeHeader(fout, fin.readline())
        else:
            fin.readline()
        count = count + 1

        for line in fin.readlines():
            fout.write(sampleName +'\t')
            data = string.split(line,'\t')
            data[1] = str(int(float(data[1])))  #start
            data[2] = str(int(float(data[1])))   #end
            fout.write(string.join(data[0:3],'\t') +'\t')
            fout.write(data[dataCol]+'\t')
            fout.write(string.join(data[3:],'\t'))
        fin.close()
    fout.close()


if len(sys.argv[:])!= 3:
    print "python parseCNV.py dir output"
    print
    sys.exit()

inDir = sys.argv[1]
if inDir[-1]!= "/":
    inDir = inDir +"/"
output = sys.argv[2]
process(inDir, output)
