import string, os, sys,stat
import json

sys.path.insert(0,"../CGDataNew")

from CGDataLib import *
import TCGAUtil

def cohort (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
        
    variable ="_cohort"
    value = "TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"
    doDerived=0
    cohort_variable (variable, value, inDir, outDir, cancer, REALRUN, doDerived)

def primary_site (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    if cancer in ["COADREAD","LUNG","PANCAN","PANCAN12","GBMLGG"]:
        return
        
    variable ="_primary_site"
    value = string.join(TCGAUtil.anatomical_origin[cancer],",")
    doDerived =1
    cohort_variable (variable, value, inDir, outDir, cancer, REALRUN, doDerived)

def primary_disease (inDir, outDir, cancer, flog,REALRUN):
    print cancer, sys._getframe().f_code.co_name
    if cancer in ["COADREAD","LUNG","PANCAN","PANCAN12","GBMLGG"]:
        return
        
    variable ="_primary_disease"
    value = TCGAUtil.cancerGroupTitle[cancer]
    doDerived =1
    cohort_variable (variable, value, inDir, outDir, cancer, REALRUN)

def cohort_variable (var, value, inDir, outDir, cancer, REALRUN, doDerived):
    print inDir
    print outDir

    if REALRUN:
        ignore =1
        bookDic=cgWalk(inDir,ignore)
        
        existMaps = collectSampleMaps(bookDic)
        missingMaps=  collectMissingSampleMaps(bookDic)

        #removeExistMaps
        for map in existMaps:
            if map not in missingMaps:
                missingMaps[map]=existMaps[map]
        
        # all aliquote uuid dic
        aliquote_dic =TCGAUtil.uuid_Aliquot_all()
        sample_dic =TCGAUtil.uuid_Sample_all()

        if len(missingMaps)!=1:
            return

        map = missingMaps.keys()[0]
        print map
        samples =[]
        for name in missingMaps[map]:
            obj=bookDic[name]
            
            if obj['type']=="genomicMatrix":
                fin =open(obj['path'],'U')
                for sample in string.split(fin.readline()[:-1],"\t")[1:]:
                    if sample =="":
                        print name, "has bad empty sample id"
                        sys.exit()
                    if sample not in samples:
                        samples.append(sample)
                fin.close()
            
            #take too long
                """
            if obj['type']=="mutationVector":
                fin =open(obj['path'],'U')
                fin.readline()
                while 1:
                    line = fin.readline()
                    if string.strip(line) =="":
                        break
                    sample = string.split(line,'\t')[0]
                    if sample not in samples:
                        samples.append(sample)
                        print sample, obj['path']
                fin.close()
            """
        intDic={}
        for sample in samples:
            #TCGA uuid handling
            uuid =sample
            TCGAbarcode =""
            if uuid[0:4]!="TCGA": 
                if aliquote_dic.has_key(string.lower(uuid)):
                    TCGAbarcode = aliquote_dic[string.lower(uuid)]
                else:
                    TCGAbarcode =  uuid
            else:
                TCGAbarcode = sample

            intID= TCGAUtil.barcode_IntegrationId(TCGAbarcode)
            if intID == None: # ids is on patient level above integration level
                continue 
            if not intDic.has_key(intID):
                intDic[intID]=""

        outfile = outDir+cancer+"/"+ var
        fout =open(outfile,"w")
        fout.write("sample\t"+var+"\n")
        for intId in intDic:
            fout.write(intId+"\t"+ value+"\n")
        fout.close()

    #data josn
    J={}
    J["version"]= datetime.date.today().isoformat()
    J["name"]="TCGA_"+cancer+"_"+var
    J["type"]= "clinicalMatrix" 
    J["dataSubType"]="phenotype"
    J[":sampleMap"]="TCGA."+cancer+".sampleMap"
    J["cohort"]="TCGA "+TCGAUtil.cancerHumanReadable[cancer]+" ("+cancer+")"

    outfile = outDir+cancer+"/"+var
    oHandle = open(outfile +".json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()

    if doDerived:
        if cancer in ["LUAD","LUSC"]:
            derived_cancer="LUNG"
            doDerivedCancer(var, outDir, cancer, derived_cancer, REALRUN)
        if cancer in ["COAD","READ"]:
            derived_cancer="COADREAD"
            doDerivedCancer(var, outDir, cancer, derived_cancer, REALRUN)
        if cancer in ["GBM","LGG"]:
            derived_cancer="GBMLGG"
            doDerivedCancer(var, outDir, cancer, derived_cancer, REALRUN)

def doDerivedCancer( var, outDir, cancer, derived_cancer, REALRUN):
    if derived_cancer =="":
        return
    
    os.system("cp "+outDir+cancer+"/"+var +" " +outDir+derived_cancer+"/"+var+"_" +cancer)

    J={}
    J["version"]= datetime.date.today().isoformat()
    J["name"]="TCGA_"+derived_cancer+"_"+var+"_"+cancer
    J["type"]= "clinicalMatrix" 
    J[":sampleMap"]="TCGA."+derived_cancer+".sampleMap"
    oHandle = open(outDir+derived_cancer+"/"+var+"_" +cancer +".json","w")
    oHandle.write( json.dumps( J, indent=-1 ) )
    oHandle.close()

    return
