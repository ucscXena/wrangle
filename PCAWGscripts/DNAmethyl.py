import os, sys

def process (inputdir):
    for file in os.path.listdir(inputdir):
        print file

if len(sys.argv[:]) != 2:
    print "python DNAmethyl.py inputdir"
    sys.exit()

inputdir = sys.argv[1]
process (inputdir)
