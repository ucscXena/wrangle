import string, os, sys
import uuid, json

def generateFilters (input,  column, states, filePrefix, cohort, outdir):
    def f1 (x): 
        if str.isalnum (x):
            return x
        else:
            return "_"

    outfile ={}
    prefix = "filter_"
    if filePrefix :
        prefix = prefix +filePrefix+"_"

    for state in states:
        cleaned =["_"]
        for i in range(0, len(state)):
            letter = state[i]
            s = f1 (letter)
            if s =="_" and cleaned[-1]=="_":
                continue
            cleaned.append(s)
        if cleaned[-1]=="_":
            cleaned = cleaned[:-1]
        cState = string.join(cleaned[1:],'')
        fout = open(prefix+cState,'w')
        fout.write("sample\n")
        outfile[state]=fout
        #json
        fout = open(prefix+cState+".json",'w')
        J= buildJson(cohort,state)
        fout.write( json.dumps( J, indent=2 ) )
        fout.close()

    fin = open(input,'r')
    fin.readline()
    for line in fin.readlines():
        data = string.split(line[:-1],'\t')
        sample = data[0]
        state =  data[column]
        if state =="":
            continue
        fout = outfile[state]
        fout.write(sample+"\n")
    fin.close()

    for key in outfile:
        fout = outfile[key]
        fout.close()

def getStates(input):
    states=[]
    fin = open(input,'r')
    for line in fin.readlines():
        state = string.strip(line)
        if state!="":
            states.append(state)
    fin.close()
    return states

def whichColumn (input, feature):
    fin = open(input,'r')
    columns = string.split(string.strip(fin.readline()),'\t')
    for i in range(1,len(columns)):
        if columns[i] == feature:
            return i
    fin.close()

def buildJson(cohort, state):
    J={}
    J["type"]="clinicalMatrix"
    J["dataSubType"]="filter"
    J["cohort"]=cohort
    J["label"]=state
    return J

if len(sys.argv[:])< 5:
    print "python buildXenaFilter.py clinicalMatrix feature cohort(use_in_json) outputDir optional_filenamePrefix(e.g.GTEX->filter_GTEX_xxx)"
    sys.exit()

input = sys.argv[1]
feature = sys.argv[2]
column = whichColumn (input, feature)
cohort  = sys.argv[3]
outdir = sys.argv[4]

if not column:
    print column, "bad feature name"
    sys.exit()

if outdir[-1]!="/":
    outdir=outdir+"/"

filePrefix =""
if len(sys.argv[:]) ==6:
    filePrefix = sys.argv[5]

output = str(uuid.uuid4())
os.system("cut -f "+str(column+1)+" " +input +" |sed 1d | sort |uniq > "+ output)
states = getStates(output)
os.system("rm -f "+output)

generateFilters (input,  column, states, filePrefix, cohort, outdir)
