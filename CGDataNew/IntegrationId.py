import string,copy

class IntegrationId():

    def __init__ (self,name,readHandle=None):
        self.__list=[]
        self.__name=name
        
        if readHandle==None:
            return
        for line in readHandle.readlines():
            id =string.strip(line)
            self.addId(id)
        return

    def getName(self):
        return self.__name
    
    def length():
        return len(self.__list)

    def getList(self):
        return copy.deepcopy(self.__list)

    def addId(self, id):
        if id is None:
            return
        if id not in self.__list:
            self.__list.append(id)
        return
    
    def store(self,oHandle):
        for id in self.__list:
            oHandle.write(id+"\n")
