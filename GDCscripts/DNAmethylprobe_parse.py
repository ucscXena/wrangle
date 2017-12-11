import string, sys

def filter_by_position ((gene, pos)):
	try:
		if int(pos) < 200:
			return True
		else:
			return False
	except:
		return False

def parse (input, output):
	fin = open(input, 'r')
	fout = open(output, 'w')
	fout.write("#id\tgene\tchrom\tchromStart\tchromEnd\tstrand\n")

	fin.readline()
	c =0 
	while 1:
		line = fin.readline()
		if line == '':
			fin.close()
			break
		data = string.split(line[:-1],'\t')
		probe = data[0]
		genes = string.split(data[5],";")
		Position_to_TSS =  string.split(data[8],";")
		Feature_Type = data[10]
		
		"""
		# Promoter CpG hyper/hypo-methylation regulate gene expression.  
		# To qualify for probeMap: probe needs to be
		# 1. upstream to the gene to downstream 200bp. ( we cover distal, proximal, and core promoter)
		# 2. Also in CpG island annotation site.

		if Feature_Type == '.':
			continue
		good_pos = filter(filter_by_position, zip(genes, Position_to_TSS))
		if len(good_pos) == 0:
			continue
		gene_list = map(lambda x: x[0], good_pos)
		genes = gene_list
		"""

		chr = data[2]
		start = data[3]
		end = data[4]

		fout.write(probe + '\t' + string.join(set(genes),',') + '\t' + chr + '\t' + start +'\t' + end +'\t' + '.\n')
		c = c+1
	print c
	fout.close()

if len(sys.argv[:])!= 3:
	print "python DNAmethylaprobe_parse.py GDC_methyl_input probeMap\n"
	sys.exit()

parse (sys.argv[1], sys.argv[2])

