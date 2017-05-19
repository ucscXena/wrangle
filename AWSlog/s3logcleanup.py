import string, os, sys
import dateutil.parser

#unify two types of load balancer logs
#http://docs.aws.amazon.com/elasticloadbalancing/latest/application/load-balancer-access-logs.html#access-log-entry-format
#http://docs.aws.amazon.com/elasticloadbalancing/latest/classic/access-log-collection.html#access-log-entry-format
def unifyLog(mergedlog, unifiedlog):
	fin = open(mergedlog, 'r')
	fout = open(unifiedlog, 'w')
	for line in fin.readlines():
		data = string.split(line, ' ')
		if data[0] in ["http", "https", "h2", "ws", "wss"]:
			fout.write(string.join(data[1:], ' '))
		else:
			fout.write(line)
	fin.close()
	fout.close()

#remove extra lines in log file  (sorted) that are from the same ip for the same file download that are within 1 min, keep only one line
def removeVeryCloseSameIPSameFile (logfile, cleanedlog):
	# collect info
	fin = open(logfile, 'r')
	time_dic = {} #key = ip_target, values: timestamp
	
	for line in fin.readlines():
		data = string.split(line, ' ')
		try:
			timestamp = data[0]
			ip = string.split(data[2],":")[0]
			target = data[12]
		except:
			continue
		if string.find(target,'download') == -1:
			continue
		if string.find (target, '/download/meta/') != -1:
			continue
		if string.find (target, '.json') != -1:
			continue
		key = ip+"_" + target
		if key not in time_dic:
			time_dic[key]=[timestamp]
		else:
			time_dic[key].append(timestamp)
	fin.close()

	# find bad lines
	bad_line_dic = {} #key = timestamp
	for key in time_dic:
		if len(time_dic[key]) == 1:
			continue
		time_list = time_dic[key]
		for i in range (0, len(time_list) -1 ):
			d1 = dateutil.parser.parse(time_list[i])
			d2 = dateutil.parser.parse(time_list[i+1])
			d = d2 - d1
			if (d.total_seconds() < 60.0): 
				if time_list[i+1] not in bad_line_dic:
					bad_line_dic[time_list[i+1]] = [key]
				else:
					bad_line_dic[time_list[i+1]].append(key)
	
	#remove bad lines
	fin = open(logfile, 'r')
	fout = open(cleanedlog, 'w')
	for line in fin.readlines():
		data = string.split(line, ' ')
		timestamp = data[0]
		if timestamp not in bad_line_dic:
			fout.write(line)
			continue
		ip = string.split(data[2],":")[0]
		target = data[12]
		if ip+'_'+ target not in bad_line_dic[timestamp]:
			fout.write(line)

	fin.close()
	fout.close()

def process (indir):
	for root, dirs, files in os.walk(indir):
		for file in files:
			filepath = root + '/'+ file
			if filepath[-7:] == ".log.gz":
				os.system("gunzip " + filepath)

	for root, dirs, files in os.walk(indir):
		for dir in dirs:
			foundLog = 0
			for filepath in os.listdir(root + '/' + dir):
				if filepath[-4:] == ".log":
					foundLog = 1
			if foundLog == 0:
				continue

			dirpath = root+'/'+dir
			mergedlog = root+'/'+dir + '/rawMerge'  
			unifiedlog = root+'/'+dir + '/unifiedMerge'  
			sortedlog = root+'/'+dir + '/sortMerge'  
			cleanedlog = root+'/'+dir + '/cleanMerge'  
			print "cleaning", dirpath
			os.system("cat " + dirpath + "/*.log > " + mergedlog)
			unifyLog(mergedlog, unifiedlog)
			os.system("LC_ALL='C' sort " + unifiedlog + " > " + sortedlog)
			removeVeryCloseSameIPSameFile (sortedlog, cleanedlog)
			os.system("rm " + string.join([mergedlog, unifiedlog, sortedlog], ' '))

if len(sys.argv[:]) != 2 :
	print "python s2logcleanup.py directory_name"
	sys.exit()

indir = sys.argv[1]
process (indir)
