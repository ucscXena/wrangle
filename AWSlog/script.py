import string,sys

fin = open(sys.argv[1],'r')
t1=0
t2=0
for line in fin.readlines():
    if string.strip(line) == "":
        continue
    data =string.split(line,' ')
    n1  = int(data[5])
    n2 = int(data[6])
    t1 = t1 + n1
    t2 = t2 + n2
fin.close()
print "IN", (t1)/1e9,"GB"
print "OUT", t2/1e9, "GB"
print "total", (t1+t2)/1e9,"GB"
