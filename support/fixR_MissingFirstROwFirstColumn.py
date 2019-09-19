import sys

def fix (input, output):
    fin = open(input,'r')
    fout = open(output,'w')
    header = fin.readline().split('\t')
    header.insert(0,'sample')
    fout.write('\t'.join(header))
    while 1:
        line = fin.readline()
        if line == '':
            break
        fout.write(line)
    fin.close()
    fout.close()
    
if len(sys.argv[:])!= 3:
    print ("python fixR_MissingFirstROwFirstColumn.py infile outfile")
    sys.exit()

input =sys.argv[1]
output = sys.argv[2]
fix (input, output)
