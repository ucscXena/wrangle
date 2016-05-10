import string, os, sys
import uuid, json

def generateFilters (input,  column, states, cohort, outdir):
    def f1 (x): 
        if str.isalnum (x):
            return x
        else:
            return "_"

    outfile ={}
    for state in states:
        cleaned =[]
        for letter in state:
            s = f1 (letter)
            if s =="_" and cleaned[-1]=="_":
                pass
            else:
                cleaned.append(s)
        cState = string.join(cleaned,'')
        fout = open(outdir+"filter_"+cState,'w')
        fout.write("sample\n")
        outfile[state]=fout
        #json
        fout = open(outdir+"filter_"+cState+".json",'w')
        J= buildJson(cohort,state)
        fout.write( json.dumps( J, indent=-1 ) )
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

if len(sys.argv[:])!= 5:
    print "python buildXenaFilter.py clinicalMatrix feature cohort(use_in_json) outputDir"
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

output = str(uuid.uuid4())
os.system("cut -f "+str(column+1)+" " +input +" |sed 1d | sort |uniq > "+ output)
states = getStates(output)
os.system("rm -f "+output)

generateFilters (input,  column, states, cohort, outdir)
