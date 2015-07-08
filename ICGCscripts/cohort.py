import string,os,sys
import icgcLib

def process(info, infile, fout):
    project = os.path.basename(infile)[9:].replace(".tsv","")
    cohort=  info[project]["xenaSuffix"]
    primary_disease = info[project]["primary_disease"]

    fin= open(infile,'rU')
    fin.readline()
    for line in fin.readlines():
        if line =="":
            break
        sample = string.split(line,"\t")[0]
        fout.write(sample+"\t"+cohort+"\t"+primary_disease+"\n")
    return

if __name__ == '__main__' :
    if len(sys.argv[:]) <4:
        print "python mergeClinicalMatrixFiles.py output inputfile(s)"
        print "this is merging data A+B=C\n"
        sys.exit()

    inFiles = sys.argv[2:]
    print inFiles
    outfile = sys.argv[1]
    print outfile
    
    info =icgcLib.getCohortInfo() #xenaSuffix primary_disease

    fout = open(outfile,"w")
    fout.write("Sample\t_Cohort\t_primary_disease\n")

    for infile in inFiles:
        process(info, infile, fout)
    fout.close()
