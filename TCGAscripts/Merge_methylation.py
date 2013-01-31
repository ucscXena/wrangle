import string, os, sys
import json,datetime

LEVEL="data.Level_3"
BETAOFFSET =-0.5

import TCGAUtil
sys.path.insert(0,"../CGDataNew")
from CGDataUtil import *

#gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2011102600.0.0.tar.gz

def Merge_methylation27(inDir, outDir, cancer,flog,REALRUN):
    PATHPATTERN= "Merge_methylation__humanmethylation27"
    Merge_methylation(inDir, outDir, cancer,flog, PATHPATTERN,REALRUN)
    return
        
def Merge_methylation(inDir, outDir, cancer,flog,PATHPATTERN,REALRUN):
    if cancer not in ["COADREAD"]:
        return
    print cancer, __name__

    garbage=["tmptmp/"]
    if os.path.exists( "tmptmp/" ):
        os.system("rm -rf tmptmp/*")
    else:
        os.system("mkdir tmptmp/")

    #figure out the FH date
    if inDir[-1]!="/":
        s = inDir+"/"
    else:
        s=inDir
    FHdate = string.split(s,"/")[-2]

    #single file in dir mode
    dataDir =""
    lastDate=None
    for file in os.listdir(inDir):
        #find the file
        if string.find(file,PATHPATTERN)!=-1 and string.find(file,LEVEL)!=-1 and string.find(file,"md5")==-1:
            pass
        else:
            continue
        
        if not os.path.exists(inDir +file+".md5"):
            print "file has no matching .md5 throw out", file
            continue

        #file date
        lastDate=  datetime.date.fromtimestamp(os.stat(inDir+file).st_mtime)


        #is tar.gz?, uncompress
        if string.find(file,".tar.gz")!=-1:
            os.system("tar -xzf "+inDir+file +" -C tmptmp/") 
            dataDir ="tmptmp/"+string.replace(file,".tar.gz","")+"/"
            break

    #make sure there is data
    if lastDate == None:
        cleanGarbage(garbage)
        print "ERROR expect data, but no data in dirpath", dataDir, cancer, __name__
        return

    if REALRUN and (dataDir =="" or not os.path.exists(dataDir)):
        cleanGarbage(garbage)
        print "ERROR expect data, but wrong dirpath", dataDir, cancer, __name__
        return


    #set output dir
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if not os.path.exists( outDir +cancer+"/"):
        os.makedirs( outDir+cancer+"/" )

    #data processing
    for file in os.listdir(dataDir):
        print file
        pattern ="data"
        if string.find(file,pattern)!=-1:
            cgFileName= string.split(dataDir[7:],PATHPATTERN)[0]+PATHPATTERN
            cgFileName = cgFileName
            infile = dataDir+file

            if REALRUN:
                outfile = outDir+cancer+"/"+cgFileName
                average = methylationProcess(infile,outfile, BETAOFFSET)
                print average, BETAOFFSET
    
            if not REALRUN:
                iHandle = open(outDir+cancer+"/"+cgFileName+".json","r")
                J = json.loads(iHandle.read())
                iHandle.close()
                average = J["mean"]

            oHandle = open(outDir+cancer+"/"+cgFileName+".json","w")

            J={}
            #stable
            if PATHPATTERN =="Merge_methylation__humanmethylation27":
                namesuffix="HumanMethylation27_FH"
                suffix="HumanMethylation27 Firehose"
                J["PLATFORM"]= "Illumina Infinium HumanMethylation27"
                J["shortTitle"]="DNA Methylation (Methylation27_Firehose)"

            J["cgDataVersion"]=1
            J["longTitle"]="TCGA "+TCGAUtil.cancerOfficial[cancer]+" ("+cancer+") DNA methylation ("+suffix+")"
            J[":dataSubType"]="DNAMethylation"
            J["redistribution"]= True
            J["groupTitle"]="TCGA "+TCGAUtil.cancerGroupTitle[cancer]
            #J["priority"]= TCGAUtil.browserPriority[cancer]
            J["dataProducer"]= "TCGA FIREHOSE pipeline"
            J["url"]= "http://gdac.broadinstitute.org/runs/stddata__"+FHdate[0:4]+"_"+FHdate[4:6]+"_"+FHdate[6:8]+"/data/"+cancer+"/"+FHdate[0:8]+"/gdac.broadinstitute.org_"+cancer+"."+ PATHPATTERN+"__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3."+FHdate+".0.0.tar.gz"
            #J["url"]=TCGAUtil.remoteBase \
            #          +string.replace(inDir,TCGAUtil.localBase,"") \
            #          + string.replace(dataDir,"tmptmp/","")[:-1]+".tar.gz"
            J["version"]= datetime.date.today().isoformat()
            J["wrangler"]= "cgData TCGAscript "+ __name__ +" processed on "+ datetime.date.today().isoformat()
            J["valOffset"] = str(BETAOFFSET)
            J["mean"]= average
            #change description
            J["gain"]=2.0
            J["wrangling_procedure"]= "FIREHOSE data download from TCGA DCC, processed at UCSC into cgData repository"
            J["description"]= "The dataset shows TCGA "+ TCGAUtil.cancerOfficial[cancer]+" ("+cancer+")"\
                              " DNA methylation data."+ \
                              " DNA methylation profile was measured experimentally using the "+J["PLATFORM"]+" platform. Beta values were derived at the Johns Hopkins University and University of Southern California TCGA genome characterization center. DNA methylation values, described as beta values, are recorded for each array probe in each sample via BeadStudio software. DNA methylation beta values are continuous variables between 0 and 1, representing the ratio of the intensity of the methylated bead type to the combined locus intensity. Thus higher beta values represent higher level of DNA methylation, i.e. hypermethylation and lower beta values represent lower level of DNA methylation, i.e. hypomethylation.\n\n We observed a bimodal distribution of the beta values from both methylation27 and methylation450 platforms, with two peaks around 0.1 and 0.9 and a relatively flat valley around 0.2-0.8.  The bimodal distribution is far more pronounced and balanced in methylation450 than methylation27 platform. In the methylation27 platform, the lower beta peak is much stronger than the higher beta peak, while the two peaks are of similar hight in the methylation450 platform. During data ingestion to UCSC cgData repository, the beta values were offset by " + str(BETAOFFSET)+" to shift the whole dataset to values between -0.5 to +0.5. The average of the unshifted beta values of this dataset is "+str(average)+\
                              ", thus much of the heatmap appears hypomethylated (blue)."+\
                              " Microarray probes are mapped onto the human genome coordinates using cgData probeMap derived from GEO GPL8490 record."
            J["description"]=J["description"]+" Here is a <a href=\"http://www.illumina.com/documents/products/appnotes/appnote_dna_methylation_analysis_infinium.pdf\" target=\"_blank\"><u>reference</u></a> to Illumina Infinium BeadChip DNA methylation platform beta value."

            J["description"] = J["description"] +"<br><br>"+TCGAUtil.clinDataDesc
                
            #change cgData
            J["name"]="TCGA_"+cancer+"_"+namesuffix
            name = trackName_fix(J['name'])
            if name ==False:
                message = "bad object name, need fix otherwise break database after loading, too long "+J["name"]
                print message
                flog.write(message+"\n")
                return
            else:
                J["name"]=name        
            if PATHPATTERN =="Merge_methylation__humanmethylation27":
                J[":probeMap"]= "illuminaMethyl27K_gpl8490"
            J["type"]= "genomicMatrix" 
            J[":sampleMap"]="TCGA."+cancer+".sampleMap"
            oHandle.write( json.dumps( J, indent=-1 ) )
            oHandle.close()
            
    cleanGarbage(garbage)
    return

def cleanGarbage(garbageDirs):
    for dir in garbageDirs:
        os.system("rm -rf dir")
    return

def methylationProcess(infile, outfile, offset):
    ID =0
    BETASTART =1
    SEG =4

    total=0.0
    count =0
    fin = open(infile,"r")
    fin.readline()
    fin.readline()
    for line in fin.readlines():
        data = string.split(string.strip(line[:-1]),"\t")
        for i in range (BETASTART, len(data), SEG):
            value = data[i]
            if value in ["","na","NA","null","NULL"]:
                continue
            value = float(value)
            total = total +value
            count = count +1
    average = total /count
    fin.close()

    fin = open(infile,"r")
    fout= open(outfile,"w")
    #header
    line= fin.readline()
    data = string.split(string.strip(line[:-1]),"\t")
    fout.write(data[ID])
    for i in range(BETASTART, len(data), SEG):
        fout.write("\t"+data[i])
    fout.write("\n")

    fin.readline()
    #data
    for line in fin.readlines():
        data = string.split(string.strip(line[:-1]),"\t")
        fout.write(data[ID])
        for i in range(BETASTART, len(data), SEG):
            value = data[i]
            if value in ["","na","NA","null","NULL"]:
                fout.write("\tNA")
            else:
                value = float(value)
                value = value +offset
                fout.write("\t"+"%.4f" % (value))
        fout.write("\n")
    fin.close()
    fout.close()
    return average
