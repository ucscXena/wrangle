import string
import copy
import re
import csv

from CGDataUtil import *
import ClinicalFeatureNew

class ClinicalMatrixNew():
    def __init__ (self, rFHandle,name,FirstColAuto=False, ClinF=None, SkipLines=[], AllowDupCOL = False):
        #1stColAuto== True, replace first column with autoincrement index 1,2,3,...
        # return emptySelf if fail to initiate
        self.__name=""
        self.__ROW=0
        self.__COL=0
        self.__COLs=[] # col header
        self.__ROWs=[] # row header
        self.__DATA={}
        self.__maxColLen =64  # maximum length of column name
        self.__IntToCategoryThreshold =10 # for int features, if <=max -> category, if > max ->float
        emptySelf =copy.deepcopy(self)
        self.__name=name

        if rFHandle== None:
            return
        try:
            rFHandle+"a"
            readHandle=open(rFHandle,'r')
        except TypeError:
            readHandle=rFHandle

        # load header validate header COLs <self.__maxColLen characters
        lineReader= csv.reader(readHandle,delimiter='\t', quotechar='"')


        while 1:
            data = lineReader.next()
            if lineReader.line_num in SkipLines:
                continue
            else:
                break

        ignoreCol=[]
        for i in range (1,len(data)):
            d= data[i]
            d =string.strip(d)
            if d=="":
                ignoreCol.append(i)
                continue

            #fix col
            new_d= col_fix(d)
            if new_d!=d:
                #print "WARNING, feature name is modified", d, "new name:", new_d
                # modify ClinF if exist
                if ClinF:
                    #ClinF.replaceFeatureName(d, new_d)
                    ClinF.replicateFeatureName(d, new_d)
                d =new_d
                
            #do not use _PATIENT
            """
            if d in ["_PATIENT","_INTEGRATION"]:
                ignoreCol.append(i)
                continue
            """
            if d not in self.__COLs:
                self.__COLs.append(d)
            else:
                if AllowDupCOL:
                    ignoreCol.append(i)
                else:
                    return
                
        self.__COL = len(self.__COLs)


        # load data
        for line in lineReader:
            if lineReader.line_num in SkipLines:
                continue

            data =line
            #bad line
            if ignoreCol==[] and len(data[1:]) != self.__COL:
                print "WARNING detected bad line sample=",data[0]
                print len(data[1:]), self.__COL, len(ignoreCol)
                for i in range (0, self.__COL-len(data[1:])):
                    data.append("")
                print "Fix", len(data[1:]), self.__COL, len(ignoreCol)

            if ignoreCol!=[]  and len(data[1:]) != self.__COL+len(ignoreCol):
                print "WARNING detected bad line sample=",data[0]
                print len(data[1:]), self.__COL, len(ignoreCol)
                for i in range (0, self.__COL+len(ignoreCol) -len(data[1:])):
                    data.append("")
                print "Fix", len(data[1:]), self.__COL+len(ignoreCol), len(ignoreCol)
                continue
            
            for i in range(0,len(data)):
                data[i]= string.strip(data[i])
            if FirstColAuto:
                sample=str(lineReader.line_num)
            else:
                sample=data[0]

            if sample not in self.__ROWs:
                self.__ROWs.append(sample)
                self.__ROW = self.__ROW+1
            if not self.__DATA.has_key(sample):
                self.__DATA[sample]={}
            else:
                print "WARNING detected duplicate sample", sample, rFHandle,name
                continue
            
            countIgnore=0

            for i in range (1, len(data)):
                if i in ignoreCol:
                    countIgnore=countIgnore+1
                    continue
                elif ignoreCol!=[]:
                    col = self.__COLs[i-1-countIgnore]
                else:
                    col = self.__COLs[i-1]

                if not self.__DATA[sample].has_key(col) :
                    if data[i] in ["NA","na","NULL", "null"]: #bad formatting in clinical data, clinical data NULL use empty
                        data[i]=""
                    self.__DATA[sample][col]=data[i]
                else:
                    self = emptySelf
                    return  

        if self.validate() != True:
            self = emptySelf
            return

        readHandle.close()


    def getName(self):
        return self.__name

    def getROWnum(self):
        return self.__ROW

    def getCOLnum(self):
        return self.__COL

    def getROW(self, row):
        if row not in self.__ROWs:
            return None
        return copy.deepcopy(self.__DATA[row])

    def getDATA(self, row, col):
        if row not in self.__ROWs:
            return None
        if col not in self.__COLs:
            return None
        return self.__DATA[row][col]

    def setDATA(self, row, col, value):
        if row not in self.__ROWs:
            return False
        if col not in self.__COLs:
            return False
        self.__DATA[row][col] = value
        return True

    def getROWs(self):
        return self.__ROWs[:]

    def hasRow(self, row):
        if row in self.__ROWs:
            return True
        else:
            return False
        
    def getCOLs(self):
        return self.__COLs[:]

    def getColStates (self, col):
        if col not in self.__COLs:
            return None
        states =[]
        for row in self.__ROWs:
            val = self.__DATA[row][col]
            if val!="" and val not in states:
                states.append(val)
        return states

    def removeCols(self, cols, validation=False):
        for col in cols:
            if col not in self.__COLs:
                continue
            self.__COL=self.__COL-1
            self.__COLs.remove(col)

            for row in self.__ROWs:
                self.__DATA[row].pop(col)
        if validation:
            return self.validate()
        else:
            return True

    def onlyKeepRows(self, rows, validation=False):
        oldRows = self.__ROWs[:]
        for row in oldRows: 
            if row not in rows:
                self.__ROW = self.__ROW-1
                self.__ROWs.remove(row)
                self.__DATA.pop(row)
        if validation:
            return self.validate()
        else:
            return True

    def removeRows(self, rows, validation=False):
        oldRows = self.__ROWs[:]
        for row in rows: 
            if row in oldRows:
                self.__ROW = self.__ROW-1
                self.__ROWs.remove(row)
                self.__DATA.pop(row)
        if validation:
            return self.validate()
        else:
            return True

    #function to add a new col as the last column
    def addOneColWithSameValue(self, colName, value):
        if colName in self.__COLs:
            return False
        
        #build out uniq new cols
        self.__COL = self.__COL+ 1
        self.__COLs.append(colName)
        for sample in self.__ROWs:
            self.__DATA[sample][colName]=value

        return True

    #function to add a new col from an existing col as the last column
    def addNewColfromOld(self, newCol, oldCol):
        if oldCol not in self.__COLs:
            return False
        if newCol in self.__COLs:
            return False
        
        #build out uniq new cols
        self.__COL = self.__COL+ 1
        self.__COLs.append(newCol)
        for sample in self.__ROWs:
            self.__DATA[sample][newCol]= self.__DATA[sample][oldCol]
        return True

    def addNewCols(self,aCMatrix,validation=False):
        print self.__name, self.__ROW
        print aCMatrix.getName(), aCMatrix.getROWnum()

        uniqCols=[]
        dupCols=[]
        for col in aCMatrix.getCOLs():
            #to do the longer version should go to the clinicalFeature file
            if len(col) >self.__maxColLen:
                newCol = col_fix(col)
                print "WARNING, feature name >",self.__maxColLen, " characters", col, newCol
                col=newCol

            states =aCMatrix.getColStates (col)
            type = aCMatrix.isTypeCategory (col)
            if type in ["category"]:
                if col in self.__COLs and colStringLength(states+self.getColStates(col)) > self.__maxStateTotalLen:
                    print "WARNING feature states length > "+ str(self.__maxStateTotalLen)+" throw out", col, self.__name
                    continue
                elif colStringLength(states) > self.__maxStateTotalLen:
                    print "WARNING feature states length > "+ str(self.__maxStateTotalLen)+" throw out", col, self.__name
                    continue

            if col not in self.__COLs:
                uniqCols.append(col)
            else:
                dupCols.append(col)
        origSelf= copy.deepcopy(self)

                        
        #build out uniq new cols
        self.__COL = self.__COL+ len(uniqCols)
        self.__COLs.extend(uniqCols)
        tempDic={}
        for col in aCMatrix.getCOLs():
            if col in uniqCols:
                tempDic[col]=""
        for sample in self.__ROWs:
            self.__DATA[sample].update(tempDic)

        #build out new sample rows  if needed
        missingSample=[]
        for sample in aCMatrix.getROWs():
            if not self.__DATA.has_key(sample):
                missingSample.append(sample)
        self.__ROW = self.__ROW +len(missingSample)
        self.__ROWs.extend(missingSample)

        tempDic={}        
        for col in self.__COLs:
            tempDic[col]=""
        for sample in missingSample:
            self.__DATA[sample]= copy.deepcopy(tempDic)
            
        for sample in aCMatrix.getROWs():
            for col in dupCols:
                if aCMatrix.getDATA(sample,col) is None or  string.strip(aCMatrix.getDATA(sample,col))=="":
                    continue
                if self.__DATA[sample][col]!=aCMatrix.getDATA(sample,col):
                    if self.__DATA[sample][col] is None or string.strip(self.__DATA[sample][col])=="":
                        #if col=="days_to_new_tumor_event_after_initial_treatment":
                        #    print sample,col, self.__DATA[sample][col], string.strip(aCMatrix.getDATA(sample,col)), "choose", string.strip(aCMatrix.getDATA(sample,col))
                        self.__DATA[sample][col]=string.strip(aCMatrix.getDATA(sample,col))
                    else:
                        #if col=="days_to_new_tumor_event_after_initial_treatment":
                        #    print "WARNING use old matrix",sample,col,self.__name,"old:",self.__DATA[sample][col] ,"new:",aCMatrix.getDATA(sample,col), "choose:", self.__DATA[sample][col]
                        pass
                        
        # load uniq col aCMatrix data to self
        for sample in aCMatrix.getROWs():
            for col in uniqCols:
                self.__DATA[sample][col]= aCMatrix.getDATA(sample,col)

        if validation:
            if self.validate()!= True:
                self = origSelf
                return False
            else:
                return True
        else:
            return True
                
    # add more __ROW, but keep __COL the same
    def addNewRows(self, idlist, dataDic,validation=False):
        newIds=[]
        for id in idlist:
            if id =="" or id in self.__ROWs:
                continue
            else:
                newIds.append(id)

        testKeys = dataDic.keys()
        #testKeys.sort()
        #self.__COLs.sort()
        #if testKeys != self.__COLs:
        #    return False
        
        for col in testKeys:
            if col not in self.__COLs:
                return False
        for col in self.__COLs:
            if col not in testKeys:
                dataDic[col]=""
            
        self.__ROWs.extend(newIds)
        self.__ROW = self.__ROW+len(newIds)
        for id in newIds:
            self.__DATA[id]=copy.deepcopy(dataDic)

        if not validation:
            return True
        else:
            return self.validate()
        
    def addColRoot(self, sampleMap):
        ColRoot ="_PATIENT"
        self.removeCols([ColRoot])
        if ColRoot not in self.__COLs:
            self.__COL= self.__COL+1
            self.__COLs.append(ColRoot)

        for sample in self.__ROWs:
            self.__DATA[sample][ColRoot]=sampleMap.getRoot(sample)
        return True

    def addColIntegration(self, sampleMap,integrationList):
        ColInt ="_INTEGRATION"
        self.removeCols([ColInt])
        if ColInt not in self.__COLs:
            self.__COL= self.__COL+1
            self.__COLs.append(ColInt)
            
        for sample in self.__ROWs:
            self.__DATA[sample][ColInt]=sampleMap.getIntegrationId(sample, integrationList)
        
    """
    def isTypeCategory (self, col):
        getRoot(sample)
        return True
    """

    def isTypeCategory (self, col):
        isFloat = True
        isInt = True
        states =self.getColStates (col)
        if states == None:
            return None
        for state in states:
            if state in [ None,""]:
                continue
            try:
                float(state)
                try:
                    int(state)
                except:
                    isInt=False
            except:
                isFloat=False
                isInt =False
                break

        if len(states) <= self.__IntToCategoryThreshold:
            return True #small number of states

        if isInt and isFloat:
            return False #stay as int -> float 
        elif not isInt and isFloat:
            return False #"float"
        else:
            return True  #"category"


    def isTypeInt (self, col):
        states =self.getColStates (col)
        if states == None:
            return None
        for state in states:
            if state in [ None,""]:
                continue
            try:
                int(state)
            except:
                isInt=False
                return False
        return True

    def isTypeFloat (self, col):
        states =self.getColStates (col)
        if states == None:
            return None
        for state in states:
            if state in [ None,""]:
                continue
            try:
                float(state)
            except:
                return False
        return True


    def pushToChildren(self, parent, sMap):
        # sMap is a object SampleMapNew
        if parent not in self.__ROWs:
            return True
        children = sMap.getChildren(parent)
        missingChildren=[]
        parentData = self.getROW(parent)
        # add parent data to existing children
        for child in children:
            if child not in self.__ROWs:
                missingChildren.append(child)
                continue
            for col in self.__COLs:
                if parentData[col] in ["",None]:
                    continue
                elif self.getDATA(child, col) in ["",None]:
                    self.__DATA[child][col]= parentData[col]
                elif self.getDATA(child, col) != parentData[col]:
                    print "WARNING use old matrix:",col, \
                          "old:", child, self.getDATA(child, col) , \
                          "new:", parent,parentData[col]
                    #return False
            
        # add missingChildren
        self.addNewRows(missingChildren,parentData)
        # initiate  children pushing data
        for child in children:
            if child == parent:
                return False
            r = self.pushToChildren(child,sMap)
            if r != True:
                return False
        return True
    
    def validate (self):
        for col in self.__COLs:
            if self.__COLs.count(col) !=1:
                return False
            if len(col)>self.__maxColLen:
                return False
                
        if self.__ROW != len(self.__ROWs):
            return False
        
        for row in self.__ROWs:
            if self.__ROWs.count(row) !=1:
                return False
        if self.__DATA.keys().sort() != self.__ROWs.sort():
            return False
    
        for sample in self.__DATA.keys():
            if self.__COL != len(self.__DATA[sample]):
                return False
            self.__COLs.sort()
            testKeys =self.__DATA[sample].keys()
            testKeys.sort()
            if self.__COLs!= testKeys:
                return False

        #test all states length for only categorical cols
        for col in self.__COLs:
            isFloat =True
            states =self.getColStates (col)
            if states == None:
                return False
            type = self.isTypeCategory (col)
            if type == None:
                return False

        return True

    def badCols(self):
        badCols=[]
        for col in self.__COLs:
            states =self.getColStates (col)
            if len(states) ==0:
                badCols.append(col)
            elif len(states)==1 and states[0]=="":
                badCols.append(col)
        self.removeCols(badCols)
        return badCols


    def findBadColsNotRemove(self):
        badCols=[]
        for col in self.__COLs:
            states =self.getColStates (col)
            if len(states) ==0:
                badCols.append(col)
            elif len(states)==1 and states[0]=="":
                badCols.append(col)
        return badCols

    #repalce value in the whole matrix
    def replaceValue(self, old, new):
        for sample in self.__ROWs:
            for feature in self.__COLs:
                if string.find(self.__DATA[sample][feature],old)!=-1:
                    self.__DATA[sample][feature]=string.replace(self.__DATA[sample][feature],old,new)
        return

    #repalce value when value=cell in the whole matrix
    def replaceValueWhole(self, old, new):
        for sample in self.__ROWs:
            for feature in self.__COLs:
                if self.__DATA[sample][feature] ==old:
                    self.__DATA[sample][feature]=new
        return

    #repalce value in a specific col
    def replaceValueInCol(self, feature, old, new):
        if feature not in self.__COLs:
            return
        for sample in self.__ROWs:
            if self.__DATA[sample][feature] == old:
                self.__DATA[sample][feature]= new
        return

    def replaceColName(self,feature, newName):
        if feature not in self.__COLs:
            return
        self.__COLs.append(newName)
        for sample in self.__ROWs:
            self.__DATA[sample][newName]=self.__DATA[sample][feature]
            self.__DATA[sample].pop(feature)
        self.__COLs.remove(feature)
        self.validate()
        return
    
    def store(self, oHandle, validation=False):
        if validation and self.validate()!=True:
            return False

        # output clinicalMatrix Data
        oHandle.write("sampleID\t"+string.join(self.__COLs,"\t")+"\n")
        for sample in self.__ROWs:
            oHandle.write(sample)
            for feature in self.__COLs:
                oHandle.write("\t"+ self.__DATA[sample][feature])
            oHandle.write("\n")
        return True

    def storeSkip1stCol(self, oHandle, validation=False):
        if validation and self.validate()!=True:
            return False
        # output clinicalMatrix Data
        oHandle.write("sampleID\t"+string.join(self.__COLs[1:],"\t")+"\n")
        for sample in self.__ROWs:
            for feature in self.__COLs[:-1]:
                oHandle.write(self.__DATA[sample][feature]+"\t")
            oHandle.write(self.__DATA[sample][self.__COLs[-1]]+"\n")
        return True
