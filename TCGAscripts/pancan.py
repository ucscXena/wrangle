import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def pancan_miRNA(dir, outDir, cancer,log, REALRUN):
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_miRNA_PANCAN"
    shortTitle="PANCAN miRNA"
    longTitle="_PANCAN AWG miRNA subtype k=15 (syn2027079)"
    maxSubtype=15
    dataProducer="https://www.synapse.org/#!Synapse:2027079"
    
    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()

    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)

def pancan_DNAmethyl (dir, outDir, cancer,log, REALRUN):
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_DNAMethyl_PANCAN"
    shortTitle="PANCAN DNA methylation"
    longTitle="_PANCAN AWG DNA methylation subtype k=19 (syn1875816)"
    maxSubtype=19
    dataProducer="https://www.synapse.org/#!Synapse:1875816"
    
    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)

def pancan_mutation (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_mutation_PANCAN"
    shortTitle="PANCAN mutation"
    longTitle="_PANCAN AWG mutation subtype k=14 (syn2492003)"
    maxSubtype=14
    dataProducer="https://www.synapse.org/#!Synapse:syn2492003"
    
    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)

def pancan_PARADIGM (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_PARADIGM_ConsensusClusters_PANCAN_K12"
    shortTitle="PANCAN PARADIGM"
    longTitle="_PANCAN PARADIGM subtype k=12 (syn1807631)"
    maxSubtype=12
    dataProducer="https://www.synapse.org/#!Synapse:synn1807631"

    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype)

def pancan_UNC_RANseq (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_UNC_RNAseq_PANCAN_K16"
    shortTitle="PANCAN RNAseq"
    longTitle="_PANCAN AWG RNAseq subtype k=16 (syn1715788)"
    maxSubtype=16
    dataProducer="https://www.synapse.org/#!Synapse:syn1715788"

    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)


def pancan_cluster_of_cluster (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_Cluster_Cluster_PANCAN"
    shortTitle="PANCAN COCA"
    longTitle="_PANCAN AWG COCA subtype assignment (syn2487022)"
    maxSubtype=13
    dataProducer="https://www.synapse.org/#!Synapse:syn2487022"

    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)

def pancan_RPPA (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_RPPA_PANCAN_K8"
    shortTitle="PANCAN protein"
    longTitle="_PANCAN AWG protein RBN subtype k=8 (syn1756922)"
    maxSubtype=8
    dataProducer="https://www.synapse.org/#!Synapse:syn1756922"
    
    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype , code)

def  pancan_CNA_named (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_CNA_PANCAN_K8"
    shortTitle="PANCAN CNA k=8"
    longTitle="_PANCAN copy number subtype k=8 (syn1712142)"
    dataProducer="https://www.synapse.org/#!Synapse:syn1712142"
    
    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()
    pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, maxSubtype,code)

def individual_DNAmethyl (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name

    featureName="_PANCAN_DNAMethyl_"+ cancer 
    shortTitle= cancer+ " DNA methylation"
    longTitle="_"+cancer+" DNA methylation subtype (syn1701558)"
    dataProducer="https://www.synapse.org/#!Synapse:syn1701558"

    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()

    individual_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, code)
    
def individual_miRNA (dir, outDir, cancer,log, REALRUN):
    #print status
    print cancer, sys._getframe().f_code.co_name
    #stupid BC OVCA type
    if cancer =="OVCA":
        cancer="OV"
    featureName="_PANCAN_mirna_"+ cancer 
    shortTitle= cancer+ " mirna"
    longTitle="_"+cancer+" microRNA subtype (syn1688309)"
    dataProducer="https://www.synapse.org/#!Synapse:syn1688309"

    filepath= dir
    code=None
    if os.path.exists(filepath+".json"):
        fin = open(filepath+".json",'r')
        J = json.loads(fin.read())
        if J.has_key("code"):
            code= J["code"]
        if J.has_key("stateOrder"):
            maxSubtype = J["stateOrder"]
        fin.close()

    individual_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, code)

   
def individual_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle, dic):
    if cancer =="PANCAN":
        return
    
    filepath=dir

    if not os.path.exists(filepath) or not os.access(filepath,os.R_OK):
        return

    outfile = featureName
    fout= open(outDir+cancer+"/"+outfile,'w')

    #header
    fout.write("sample\t"+featureName+"\n")
        
    #data
    fin= open(filepath,'r')
    fin.readline()

    if dic ==None:
        fout.write(fin.read())
    else:
        for line in fin.readlines():
            id,cluster= string.split(string.strip(line),"\t")
            fout.write(id+"\t"+dic[cluster]+"\n")
    
    fout.close()
    fin.close()

    #data json
    J={}
    J['name']=featureName
    J["version"]= datetime.date.today().isoformat()
    J["type"]= "clinicalMatrix"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["cohort"]="TCGA "+TCGAUtil.cancerHumanReadable[cancer]
    J['dataProducer']=dataProducer
    J[":clinicalFeature"] = J['name']+ "_clinicalFeature"
    fout= open(outDir+cancer+"/"+outfile+".json",'w')
    fout.write( json.dumps( J, indent=-1 ) )
    fout.close()

    #clinical feature data
    fout= open(outDir+cancer+"/"+outfile+"_clinicalFeature",'w')
    fout.write(featureName+"\tshortTitle\t"+shortTitle+"\n")
    fout.write(featureName+"\tlongTitle\t"+longTitle+"\n")
    fout.write(featureName+"\tvalueType\tcategory\n")


    if TCGAUtil.featurePriority.has_key(cancer):
        if TCGAUtil.featurePriority[cancer].has_key(featureName):
            priority= TCGAUtil.featurePriority[cancer][featureName]
            fout.write(featureName+"\tpriority\t"+str(priority)+"\n")
            fout.write(featureName+"\tvisibility\ton\n")
                
    fout.close()
    
    #clinical feature json
    fout= open(outDir+cancer+"/"+outfile+"_clinicalFeature.json",'w')
    JC={}
    JC['name']= J[":clinicalFeature"]
    JC["version"]= datetime.date.today().isoformat()
    JC["type"]= "clinicalFeature"
    fout.write( json.dumps( JC, indent=-1 ) )
    fout.close()

    return


def pancan_subtype (dir, outDir, cancer,log, REALRUN, featureName,dataProducer, shortTitle, longTitle,maxSubtype, dic=None):
    filepath=dir
    if not os.path.exists(filepath) or not os.access(filepath,os.R_OK):
        return

    for cancerDir in os.listdir(outDir):
        outfile = featureName
        fout= open(outDir+cancerDir+"/"+outfile,'w')
        #header
        fout.write("sample\t"+featureName+"\n")
        
        #data
        fin= open(filepath,'r')

        fin.readline()

        if dic ==None:
            fout.write(fin.read())
        else:
            for line in fin.readlines():
                id,cluster= string.split(line[:-1],"\t")
                if cluster in ["","NA"]:
                    fout.write(id+"\t\n")
                else:
                    fout.write(id+"\t"+dic[cluster]+"\n")
        fout.close()
        fin.close()

        #data json
        J={}
        J['name']=featureName+"_"+cancerDir
        J["version"]= datetime.date.today().isoformat()
        J["type"]= "clinicalMatrix"
        J[":sampleMap"]="TCGA."+cancerDir+".sampleMap"
        J["cohort"]="TCGA "+TCGAUtil.cancerHumanReadable[cancerDir]
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

        s=""
        if isinstance(maxSubtype, int):
            for i in range (1,maxSubtype):
                fout.write(featureName+"\tstate\t"+str(i)+"\n")
                s=s+"\""+str(i)+"\","
            fout.write(featureName+"\tstate\t"+str(maxSubtype)+"\n")
            s=s+"\""+str(maxSubtype)+"\""
            fout.write(featureName+"\tstateOrder\t"+s+"\n")


        elif isinstance( maxSubtype, list):
            for i in range (0,len(maxSubtype)-1):
                fout.write(featureName+"\tstate\t"+maxSubtype[i]+"\n")
                s=s+"\""+maxSubtype[i]+"\","
            fout.write(featureName+"\tstate\t"+maxSubtype[-1]+"\n")
            s=s+"\""+maxSubtype[-1]+"\""
            fout.write(featureName+"\tstateOrder\t"+s+"\n")
            
        else:
            pass

        if TCGAUtil.featurePriority.has_key(cancerDir):
            if TCGAUtil.featurePriority[cancerDir].has_key(featureName):
                priority= TCGAUtil.featurePriority[cancerDir][featureName]
                fout.write(featureName+"\tpriority\t"+str(priority)+"\n")
                fout.write(featureName+"\tvisibility\ton\n")
                
        #clinical feature json
        fout= open(outDir+cancerDir+"/"+outfile+"_clinicalFeature.json",'w')
        JC={}
        JC['name']= J[":clinicalFeature"]
        JC["version"]= datetime.date.today().isoformat()
        JC["type"]= "clinicalFeature"
        fout.write( json.dumps( JC, indent=-1 ) )
        fout.close()

    return
