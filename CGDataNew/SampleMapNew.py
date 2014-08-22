import string,copy

class SampleMapNew():

    # DATA {} stores Data[parent]=[child, child,...]
    def __init__ (self, readHandle, name):
        # return None if fail to be validated
        self.__DATA={}
        self.__name =name

        #for only an empty sampleMap
        if readHandle==None:
            self.__DATA={}
            self.__name =name
            return
        
        for line in readHandle.readlines():
            data = string.split(string.strip(line),"\t")
            if len(data)!=2:
                self.__DATA={}
                self.__name =None
                return
            
            parent, child = data
            parent = string.strip(parent)
            child = string.strip(child)
            if parent == child:
                if not self.__DATA.has_key(parent):
                    self.__DATA[parent]=[]
            else:
                if self.__DATA.has_key(parent):
                    if child not in self.__DATA[parent]:
                        self.__DATA[parent].append(child)
                else:
                    self.__DATA[parent]=[child]
                if not self.__DATA.has_key(child):
                    self.__DATA[child]=[]

        if self.validate()!= True:
            self.__DATA={}
            self.__name =None
            return
            
    def validate(self):
        #validate that each node has only one parent
        """
        allV=[]
        for values in self.__DATA.values():
            for value in values:
                if value not in allV:
                    allV.append(value)
                else:
                    return False
        """
        # validate that each child is present in the key space
        for values in self.__DATA.values():
            for child in values:
                if not self.__DATA.has_key(child):
                    return False

        # validate that key in not in values
        for key in self.__DATA.keys():
            if key in self.__DATA[key]:
                return False
        return True

    def __len__(self):
        return len(self.__DATA)

    def getNodes(self):
        return copy.deepcopy(self.__DATA.keys())
        
    def getName(self):
        return self.__name
    
    def inData(self, id):
        return self.__DATA.has_key(id)

    def getParent(self, id):
        if not self.inData(id):
            return None
        for parent in self.__DATA.keys():
            children = self.__DATA[parent]
            if id in children:
                return parent
        return ""

    
    def isRoot(self, id):
        if not self.inData(id):
            return None 
        if self.getParent(id) =="":
            return True
        else:
            return False

    def allRoots(self):
        roots=[]
        for node in self.getNodes():
            if self.isRoot(node):
                roots.append(node)
        return roots
        
    def isLeaf(self, id):
        if not self.inData(id):
            return None 
        children = self.__DATA[id]
        if children ==[]:
            return True
        else:
            return False

    def getRoot(self,id):
        if not self.inData(id):
            return None
        parent = self.getParent(id)
        if parent == None:
            return None
        if parent =="":
            return id
        else:
            return self.getRoot(parent)

    def getIntegrationId(self, id, integrationList):
        if integrationList is None or integrationList ==[]:
            return self.getRoot(id)
        if id in integrationList:
            return id
        parent = self.getParent(id)
        if parent=="" or parent ==None:
            return ""
        else:
            return self.getIntegrationId(parent, integrationList)
        
    def getChildren(self,id):
        if not self.inData(id):
            return None
        return self.__DATA[id]
    
    def getDescendants (self,id):
        if not  self.inData(id):
            return None
        desc=self.__DATA[id]
        if desc ==[]:
            return desc
        totalDesc= desc[:]
        for sample in desc:
            totalDesc.extend(self.getDescendants(sample))
        for sample in totalDesc:
            if totalDesc.count(sample)!=1:
                print sample, totalDesc.count(sample), totalDesc,"SampleMap CRITICAL ERROR"
        return totalDesc

    def addNode(self,id):
        if self.inData(id):
            return False
        self.__DATA[id]=[]
        return True

    def addLink(self,parent,child):
        if not self.inData(parent):
            self.addNode(parent)
        if not self.inData(child):
            self.addNode(child)
        if child not in self.__DATA[parent]:
            self.__DATA[parent].append(child)
        return True
        
    def store(self,oHandle):
        for parent in self.__DATA.keys():
            children = self.__DATA[parent]
            if children ==[] and self.isRoot(parent):
                if parent !="" and  parent ==  string.strip(parent):
                    oHandle.write(parent+"\t"+parent+"\n")
            else:
                for child in children:
                    if parent !="" and child !="":
                        oHandle.write(parent+"\t"+child+"\n")
                
