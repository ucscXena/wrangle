import string, sys

def removeRowGenomicMatrix_nodata  (inputFile, outputFile, nodata_values):
    fin = open(inputFile,'U')
    fout = open(outputFile,'w')

    #header
    line = fin.readline()
    fout.write(line)

    #nCOL
    nCOL = len(string.split(line[:-1],'\t'))

    #data
    while 1:
        line = fin.readline()
        if line =="":
            break
        data = string.split(line[:-1],'\t')
        if nCOL != len(data):
            print "Wrong nCOL", data[0]
            sys.exit()

        found = 0
        for i in range(1, nCOL):
            if data[i] not in nodata_values:
                found =1
                break
        if not found:
            continue

        fout.write(line)

        #new_data = data[i][:]
        #for i in range(1, nCOL):
        #    if float(data[i]) == nodata_value:
        #        new_data[i] = ""
        #fout.write(string.join(new_data,'\t'))
        #fout.write('\n')

    fin.close()
    fout.close()

if len(sys.argv[:]) < 4:
    print "python removeRowGenomicMatrix_nodata.py input output no_data(value like 0) no_data(value like NA)"
    print
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
nodata_values = sys.argv[3:]
removeRowGenomicMatrix_nodata (inputfile, outputfile, nodata_values)
