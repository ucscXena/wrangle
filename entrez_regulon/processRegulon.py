import string,sys

fRegulon = "regulonBRCA"
fEntrezMapping = "entrez_id_name"
convertedRegulon = "regulonBRCA_hugo"

def processRegulon(id_name_dic,fRegulon, convertedRegulon):
	fin = open(fRegulon, 'U')
	fout = open(convertedRegulon, 'w')
	line = fin.readline()
	fout.write(line)
	for line in fin.readlines():
		data = string.split(line[:-1],'\t')
		try:
			Regulator = id_name_dic[data[0]]
		except:
			print data[0]
			continue
		try:
			Target = id_name_dic[data[1]]
		except:
			print data[1]
			continue
		fout.write(string.join([Regulator, Target, data[2], data[3]], '\t')+'\n')
	fin.close()
	fout.close()

def entrez_id_name_mapping(fEntrezMapping):
	dic ={}
	fin = open(fEntrezMapping, 'U')
	fin.readline()
	for line in fin.readlines():
		data = string.split(line[:-1],'\t')
		ENTREZID = data[0]
		SYMBOL = data[1]
		dic[ENTREZID] = SYMBOL
	fin.close()
	return dic

id_name_dic = entrez_id_name_mapping(fEntrezMapping)
processRegulon(id_name_dic,fRegulon, convertedRegulon)
