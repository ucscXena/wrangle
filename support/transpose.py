import string, sys, os

def process (infile, outfile):
    fout=open(outfile,'w')
    fout.close()

    #header
    fin=open(infile,'rU')
    line = fin.readline()
    totalN = len(string.split(line,'\t'))
    fin.close()

    for i in range (0, totalN):
        start = i + 1 #linux cut
        command = "cut -f " + str(start) + " " + infile + " | tr '\\n' '\\t' | sed 's/.$/\\n/' >> " + outfile
        #print command
        os.system(command)

if len(sys.argv[:])!= 3:
    print "python transpose_memEfficient.py fin fout"
    sys.exit()

infile = sys.argv[1]
outfile= sys.argv[2]

process(infile, outfile)
