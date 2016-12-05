import string, sys

def list_to_percentiles(numbers):
    dic ={}
    total = 0 
    for i in range (0, len(numbers)):
        n = numbers[i]
        if n == '':
            continue
        total = total + 1
        if n not in dic:
            dic[n]=[]
        dic[n].append(i)

    num_list = dic.keys()
    num_list.sort()  #small to large
    result = ['' for i in range(len(numbers))]    
    running_t =0 
    for rank in xrange(len(num_list)):
        n = num_list[rank]
        original_index_list = dic[n]
        current_size = len(original_index_list)
        percentRank = running_t * 100.0 / (total)
        running_t = running_t + current_size
        for index in original_index_list:
            result[index] = percentRank
    return result

def process (infile, outfile):
    fin =open(infile,'r')
    fout =open(outfile,'w')
    fout.write(fin.readline())

    while 1:
        line = fin.readline()
        if line =="":
            break
        data =string.split(line,'\t')

        id = data[0]
        fout.write(id)

        for i in range (1, len(data)):
            try:
                data[i] = float(data[i])
            except:
                data[i] = ''

        result = list_to_percentiles(data[1:])

        for item in result:
            fout.write('\t'+ str(item))
        fout.write('\n')
    fin.close()

if len(sys.argv[:])!=3:
    print "python percentRankByRow.py inputMatrix outputMatrix"
    print
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

process (infile, outfile)
