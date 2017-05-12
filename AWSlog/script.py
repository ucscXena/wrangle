import string,sys

fin = open(sys.argv[1],'r')
t1=0
t2=0
while 1:
	line = fin.readline()
	if line == '':
		break
	if string.strip(line) == "":
		continue
	data =string.split(line,' ')
	try:
		n1  = int(data[9])
		n2 = int(data[10])
		t1 = t1 + n1
		t2 = t2 + n2
	except:
		continue
fin.close()
print "IN", (t1)/1e9,"GB"
print "OUT", t2/1e9, "GB"
print "total", (t1+t2)/1e9,"GB"
