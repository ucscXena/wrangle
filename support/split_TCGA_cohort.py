import sys, json, os

project_code_file = os.path.dirname(os.path.realpath(__file__)) + "/TCGA_project_code.tsv"

def parse_project_code (project_code_file):
	project_code = {}  #key: sample ID value: project_code
	project_list = []
	fin = open(project_code_file, 'r')
	fin.readline()
	for line in fin.readlines():
		data = line[:-1].split('\t')
		sample = data[0]
		project = data[1]
		project_code[sample] = project
		if project not in project_list:
			project_list.append(project)
	fin.close()
	return project_code, project_list 

def build_json(projects, outputdir, example_json):
	cohort_name = {
		"LAML" : "TCGA Acute Myeloid Leukemia (LAML)",
		"ACC": "TCGA Adrenocortical Cancer (ACC)",
		"CHOL": "TCGA Bile Duct Cancer (CHOL)",
		"BLCA": "TCGA Bladder Cancer (BLCA)",
		"BRCA": "TCGA Breast Cancer (BRCA)",
		"CESC": "TCGA Cervical Cancer (CESC)",
		"COAD": "TCGA Colon Cancer (COAD)",
		"UCEC": "TCGA Endometrioid Cancer (UCEC)",
		"ESCA": "TCGA Esophageal Cancer (ESCA)",
		"GBM": "TCGA Glioblastoma (GBM)",
		"HNSC": "TCGA Head and Neck Cancer (HNSC)",
		"KICH": "TCGA Kidney Chromophobe (KICH)",
		"KIRC": "TCGA Kidney Clear Cell Carcinoma (KIRC)",
		"KIRP": "TCGA Kidney Papillary Cell Carcinoma (KIRP)",
		"DLBC": "TCGA Large B-cell Lymphoma (DLBC)",
		"LIHC": "TCGA Liver Cancer (LIHC)",
		"LGG": "TCGA Lower Grade Glioma (LGG)",
		"LUAD": "TCGA Lung Adenocarcinoma (LUAD)",
		"LUSC": "TCGA Lung Squamous Cell Carcinoma (LUSC)",
		"SKCM": "TCGA Melanoma (SKCM)",
		"MESO": "TCGA Mesothelioma (MESO)",
		"UVM": "TCGA Ocular melanomas (UVM)",
		"OV": "TCGA Ovarian Cancer (OV)",
		"PAAD": "TCGA Pancreatic Cancer (PAAD)",
		"PCPG": "TCGA Pheochromocytoma & Paraganglioma (PCPG)",
		"PRAD": "TCGA Prostate Cancer (PRAD)",
		"READ": "TCGA Rectal Cancer (READ)",
		"SARC": "TCGA Sarcoma (SARC)",
		"STAD": "TCGA Stomach Cancer (STAD)",
		"TGCT": "TCGA Testicular Cancer (TGCT)",
		"THYM": "TCGA Thymoma (THYM)",
		"THCA": "TCGA Thyroid Cancer (THCA)",
		"UCS": "TCGA Uterine Carcinosarcoma (UCS)"
	}

	fin = open(example_json, 'r')
	J = json.loads(fin.read())
	fin.close()

	for project in projects:
		J['cohort'] = cohort_name[project]
		fout = open(outputdir + '/' + project + '_' + outputdir +'.txt.json', 'w')		
		fout.write(json.dumps(J, indent=4))
		fout.close()

def parse_matrix_data (data_file, project_code_dic, projects, outputdir):
	fin = open (data_file, 'r')

	#open output files
	project_file_dic ={}
	for project in projects:
		fout = open(outputdir + '/' + project + '_' + outputdir +'.txt', 'w')
		fout.write('sample')
		project_file_dic[project] = fout

	#header -- sample list
	pos_project_dic={}
	line = fin.readline()
	samples = line.strip().split('\t')
	for i in range (1, len(samples)):
		sample = samples[i]
		project = project_code_dic[sample]
		fout = project_file_dic[project]
		fout.write('\t' + sample)
		pos_project_dic[i] = project

	for project in projects:
		fout = project_file_dic[project]
		fout.write('\n')

	# data
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = line[:-1].split('\t')
		gene = data[0]
		for project in projects:
			fout = project_file_dic[project]
			fout.write(gene)
		for i in range (1, len(data)):
			value = data[i]
			project = pos_project_dic[i]
			fout = project_file_dic[project]
			fout.write('\t' + value)
		for project in projects:
			fout = project_file_dic[project]
			fout.write('\n')
	fin.close()

	# close output files
	for project in projects:
		fout = project_file_dic[project]
		fout.close()

def parse_table_data (data_file, project_code_dic, projects, outputdir):
	fin = open (data_file, 'r')

	#open output files and write header
	project_file_dic ={}
	header = fin.readline()
	for project in projects:
		fout = open(outputdir + '/' + project + '_' + outputdir +'.txt', 'w')
		fout.write(header)
		project_file_dic[project] = fout

	# data
	while 1:
		line = fin.readline()
		if line == '':
			break
		data = line[:-1].split('\t')
		sample = data[0]
		try:
			project = project_code_dic[sample]
			fout = project_file_dic[project]
			fout.write(line)
		except KeyError:
			print (sample)
	fin.close()

	# close output files
	for project in projects:
		fout = project_file_dic[project]
		fout.close()

if len(sys.argv[:]) != 5:
	print ("python split_cohort.py data_file example_json table/matrix outputdir")
	print
	sys.exit()


data_file = sys.argv[1]
example_json = sys.argv[2]
datatype = sys.argv[3]
outputdir = sys.argv[4]

project_code_dic, projects = parse_project_code (project_code_file)

os.system("mkdir " + outputdir)

build_json(projects, outputdir, example_json)

if datatype == "matrix":
	parse_matrix_data (data_file, project_code_dic, projects, outputdir)
elif datatype == 'table':
	parse_table_data (data_file, project_code_dic, projects, outputdir)
