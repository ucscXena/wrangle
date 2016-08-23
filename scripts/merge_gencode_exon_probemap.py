import string, sys, os

def sort_by_id (input, output ):
    os.system("sort -k 1 "+ input + " | uniq > " + output)
    return

def process (prev_id, data_list, fout):
    if len(data_list) == 0:
        return
    if prev_id =="":
        return

    if len(data_list) == 1:
        fout.write(prev_id + "\t")
        fout.write(string.join(data_list[0],'\t'))
        fout.write("\n")
        return


    starts =[]
    ends=[]
    for item in data_list:
        end1 = int(item[2])
        end2 = int(item[3])
        if end1 > end2:
            end1 = int(item[3])
            end2 = int(item[2])
        starts.append(end1)
        ends.append(end2)

    starts.sort()
    ends.sort()

    start = str(starts[0])
    end = str(ends[-1])

    new_list = data_list[0][:]
    new_list[2] = start
    new_list[3] = end
    fout.write(prev_id + "\t")
    fout.write(string.join(new_list,'\t'))
    fout.write("\n")

def merge_same_exon_id(input, output):
    fin = open(input, 'r')
    fout = open(output, 'w')

    #same header
    fout.write(fin.readline())
    prev_id =""
    data_list =[]
    while 1:
        line = fin.readline()
        line = string.strip(line)
        if line == "":
            break
        data = string.split(line,'\t')
        id = data [0]
        if id == prev_id:
            data_list.append(data[1:])
        else:
            process (prev_id, data_list, fout)
            prev_id = id
            data_list = [data[1:]]

    fin.close()
    fout.close()

if len(sys.argv[:]) != 3:
    print "python merge_gencode_exon_probemap.py inputfile(merge-multiple-record-same-exon-id) outputfile"
    print
    sys.exit()

input = sys.argv[1]
output = sys.argv[2]

tmpfile = input + "_sortuniq"

#sort_by_id (input)
merge_same_exon_id(tmpfile, output)
