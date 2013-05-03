import string, os, sys
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def pancan_DNAmethyl (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_DNAMethyl_PANCAN_C18"
    shortTitle="DNA methyl PANCAN subtype"
    longTitle="DNA methylation PANCAN subtype k=18"
    dataProducer="https://www.synapse.org/#!Synapse:syn1807207"
    
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle)

    
def pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle):
    filepath=dir
    if not os.path.exists(filepath):
        return

    file=os.path.basename(filepath)

    for cancerDir in os.listdir(outDir):
        outfile = "_PANCAN_"+file
        fout= open(outDir+cancerDir+"/"+outfile,'w')

        #header
        fout.write("sample\t"+featureName+"\n")
        
        #data
        fin= open(filepath,'r')
        fin.readline()
        fout.write(fin.read())
        fout.close()
        fin.close()

        #data json
        J={}
        J['name']=file
        J["version"]= datetime.date.today().isoformat()
        J["type"]= "clinicalMatrix"
        J[":sampleMap"]="TCGA."+cancer+".sampleMap"
        J['dataProducer']=dataProducer
        J[":clinicalFeature"] = J['name']+ "_clinicalFeature"
        fout= open(outDir+cancerDir+"/"+outfile+".json",'w')
        fout.write( json.dumps( J, indent=-1 ) )
        fout.close()

        #clinical feature data
        fout= open(outDir+cancerDir+"/"+outfile+"_clinicalFeature",'w')
        fout.write(featureName+"\tshortTitle\t"+shortTitle+"\n")
        fout.write(featureName+"\tlongTitle\t"+longTitle+"\n")
        fout.write(featureName+"\tvalueType\tcategory\n")
        fout.close()
                   
        #clinical feature json
        fout= open(outDir+cancerDir+"/"+outfile+"_clinicalFeature.json",'w')
        JC={}
        JC['name']= J[":clinicalFeature"]
        JC["version"]= datetime.date.today().isoformat()
        JC["type"]= "clinicalFeature"
        fout.write( json.dumps( JC, indent=-1 ) )
        fout.close()

    return
