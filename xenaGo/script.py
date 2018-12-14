import string, sys, os
import json

defaultDatafile = "defaultDatasetForGeneset.json"
step2datadir = "cnv_mutation/"

def stepDir():
	for cohort in defaults:
		os.system("mkdir " + cohort)

def stepDownload():
	for cohort in defaults:
		print cohort
		# download
		host = defaults[cohort]["simple somatic mutation"]["host"]
		dataset = defaults[cohort]["simple somatic mutation"]["dataset"]
		download = host + "/download/" + dataset + ".gz"
		os.system("curl " + download + " -o " + cohort+ "/" + os.path.basename(download))
		
		host = defaults[cohort]["copy number for pathway view"]["host"]
		dataset = defaults[cohort]["copy number for pathway view"]["dataset"]
		download = host + "/download/" + dataset + ".gz"
		os.system("curl " + download + " -o " + cohort+ "/" + os.path.basename(download))
		
		# unzip
		os.system("gunzip -f " + cohort + "/*.gz")

def stepRunExpObs():
	for cohort in defaults:
		#if cohort != "CCLE_Breast":
		#	continue

		print cohort
		# expected and observed runs
		os.system("python2.7 expectedHypergeometric.py tgac.js " + cohort + " " + cohort + "/* &")
		os.system("python2.7 observed.py tgac.js " + cohort + " " + cohort + "/* &")

def stepRunChiSquare():
	for cohort in defaults:
		#if cohort != "CCLE_Breast":
		#	continue

		print cohort
		os.system("python2.7 chiSquare.py tgac.js " + step2datadir + cohort + "_pathway_expected " 
			+ step2datadir + cohort + "_pathway_observed " + cohort + "_output")
		

def stepCompileOutput():
	label ={
		4: 'Expected',
		5: 'Observed',
		7: '%_samples',
		8: '%_samples_perGene',
		9: 'log2_ObsDivExp',
		10: 'chi_squre_value', 
		11: '1-prob',
		12: 'prob', 
		13: '-log10_prob'
	}
	seg = 13
	
	for start_col in label.keys():
		print label[start_col]
		output = step2datadir[:-1] + "_" + label[start_col]
		command ='paste '
		for cohort in defaults:
			command = command + step2datadir + cohort+"_output "
		command = command + " | cut -f 1"


		for i in range (0, len (defaults)):
			command = command + ','+str(start_col + i*seg)
		command = command + ' > tmp'
		os.system(command)

		changeHeader(label[start_col], "tmp",  output)
		
		# convert to genomicMatrix orientatin
		os.system("python2.7 ~/git/wrangle/support/transpose.py \"" + output +"\" tmp")
		os.system("mv tmp " + output)

def changeHeader(label, infile,  outfile):
	fin = open(infile,'U')
	fout = open(outfile,'w')
	data = string.split(fin.readline(),'\t')
	for i in range (1, len(data)):
		data[i] = string.replace(data[i], "_"+label, '')
	fout.write(string.join(data,'\t'))
	fout.write(fin.read())
	fin.close()
	fout.close()


fdefaultdata = open(defaultDatafile, 'U')
defaults = json.loads(fdefaultdata.read())
fdefaultdata.close()

#stepDir()
#stepDownload()
#stepRunExpObs()
#stepRunChiSquare()
stepCompileOutput()