
import os
import re
import json
import functools
from zipfile import ZipFile
import sys
import hashlib
"""
CGData object style:

Every file type documented in the CGData specification has an equivilent object
to parse and manipulate the contents of that file type. For <dataType> there
should be a CGData.<dataType> object with a <CGData> class. These classes
should extend the baseObject class. For loading they implement the 'read'
function which will parse the contents of a file from a passed file handle.
"""


OBJECT_MAP = {
    'genomicSegment': ('CGData.GenomicSegment', 'GenomicSegment'),
    'genomicMatrix': ('CGData.GenomicMatrix', 'GenomicMatrix'),
    'probeMap': ('CGData.ProbeMap', 'ProbeMap'),
    'probeLoc': ('CGData.ProbeLoc', 'ProbeLoc'),
    'aliasMap' : ('CGData.AliasMap', 'AliasMap'),
    'idDAG': ('CGData.IDDag', 'IDDag'),
    'clinicalMatrix': ('CGData.ClinicalMatrix', 'ClinicalMatrix'),
    'dataSubType': ('CGData.DataSubType', 'DataSubType'),
    'assembly': ('CGData.Assembly', 'Assembly'),
    'featureDescription': ('CGData.FeatureDescription', 'FeatureDescription'),
    'refGene' : ('CGData.RefGene', 'RefGene'),
    'idList' : ('CGData.IDList', 'IDList')
}

class FormatException(Exception):

    def __init__(self, str):
        Exception.__init__(self, str)


def has_type(type_str):
    return type_str in OBJECT_MAP

def get_type(type_str):
    mod_name, cls_name = OBJECT_MAP[type_str]
    module = __import__(mod_name, globals(), locals(), [ cls_name ])
    cls = getattr(module, cls_name)
    return cls

class UnimplementedException(Exception):
    def __init__(self, str="Method not implemented"):
        Exception.__init__(self, str)

class CGObjectBase(dict):
    """
    This is the base object for CGData loadable objects.
    The methods covered in the base case cover usage meta-information
    loading/unloading and manipulation as well as zip (cgz) file access.
    """
    __format__ = None
    def __init__(self):
        self.path = None
        self.zip = None
        self.light_mode = False
        self.loaded = False
        if 'cgformat' not in self and self.__format__ is not None:
            self['cgformat'] = self.__format__
        super(CGObjectBase,self).__init__()

    def load(self, path=None, **kw):
        """
        Load a data object in from path
        """
        if path is None and self.path is not None:
            path = self.path
        if path is None:
            raise OSError( "Path not defined" ) 
        
        if self.zip is None:
            if os.path.exists(path):
                dhandle = open(path,'rU')
                self.read(dhandle, **kw)
                dhandle.close()
        else:
            z = ZipFile(self.zip)
            dhandle = z.open(self.path, 'rU')
            self.read(dhandle, **kw)
            dhandle.close()
            z.close()
            
        self.path = path
        if (os.path.exists(path + ".json")):
            mhandle = open(path + ".json",'rU')
            meta = json.loads(mhandle.read())
            meta = dict((k, v) for k, v in meta.iteritems() if v != None)
            self.update(meta)
            mhandle.close()
        self.loaded = True

    def unload(self):
        """Call to start freeing up memory"""
        self.free()
        self.loaded = False
    
    def store(self, path=None):
        """
        Store an object onto the path provided.
        Will write a path and a path.json file.
        """
        if path is None and self.path is not None:
            path = self.path
        if path is None:
            raise OSError( "Path not defined" ) 
        meta = {}
        meta.update(self)
        if 'cgformat' in meta:
            del meta['cgformat']
        mHandle = open(path + ".json", "w")
        mHandle.write(json.dumps(meta))
        mHandle.close()
        if not self.light_mode:
            self.path = path
            dhandle = open(path, "w")
            self.write(dhandle)
            dhandle.close()            
    
    def load_keyset(self, key_predicate):
        if self.path is not None:
            if self.zip is None:
                if os.path.exists(self.path):
                    dhandle = open(self.path, 'rU')
                    out = self.read_keyset(dhandle, key_predicate)
                    for a in out:
                        yield a
                    dhandle.close()
            else:
                z = ZipFile(self.zip)
                dhandle = z.open(self.path, 'rU')
                out = self.read_keyset(dhandle, key_predicate)
                for a in out:
                    yield a
                dhandle.close()
                z.close()
        
    def read_keyset(self, handle, key_predicate=None):
        raise UnimplementedException()
    
    def read(self, handle):
        """
        The read method is implemented by the subclass that 
        inherits from CGObjectBase. It is passed a handle 
        to a file (which may be on file, in a compressed object, or
        from a network source). The implementing class then uses his handle
        to populate it's data structures. 
        """
        raise UnimplementedException()
    
    def write(self, handle):
        """
        The write method is implemented by the subclass that 
        inherits from CGObjectBase. It is passed a handle to an 
        output file, which it can use 'write' method calls to emit
        it's data.
        """
        raise UnimplementedException()

    def get_name(self):
        """
        Get object name
        """
        return self.get( 'cgdata', {} ).get( 'name', None )
    
    def get_type(self):
        """
        Get object type
        """
        return self.get('cgdata', {}).get('type', None)
    
    def get_link_map(self):
        """
        Get a dict that represents the declared file relationships from the meta-info
        """
        out = {}
        if "cgformat" in self:
            if "links" in self["cgformat"]:
                for field in self['cgformat']['links']:
                    if field in self['cgdata']:
                        if isinstance(self['cgdata'][field], str) or isinstance(self['cgdata'][field], unicode) :
                            out[field] = { 'type' : field, 'name' : self['cgdata'][field] }
                        else:
                            out[field] = { 'type' : self['cgdata'][field]['type'], 'name' : self['cgdata'][field]['name'] }
                            
        for e in ['columnKeySrc', 'rowKeySrc' ]:
            if e in self['cgdata']:
                if e not in out:
                    out[e] = {}
                link = self['cgdata'][e]
                out[e] = { 'type' : link['type'], 'name' : link['name'] }
        return out

    def add_history(self, desc):
        if not 'history' in self:
            self[ 'history' ] = []
        self[ 'history' ].append( desc )


class CGDataMatrixObject(CGObjectBase):
        
    def __init__(self):
        CGObjectBase.__init__(self)
    
    
    def get_col_namespace(self):
        """
        Return the name of the column namespace
        """
        raise UnimplementedException()

    def get_row_namespace(self):
        """
        Return the name of the row namespace
        """
        raise UnimplementedException()
    
    def get_col_list(self):
        """
        Returns names of columns
        """
        raise UnimplementedException()
        
    def get_row_list(self):
        """
        Returns names of rows
        """
        raise UnimplementedException()
    
    def get_row_map(self):
        """
        Returns map of row name indexes
        """
        raise UnimplementedException()
         
    def get_col_map(self):
        """
        Returns map of row name indexes
        """
        raise UnimplementedException()
         
    
    def get_row_pos(self, row):
        raise UnimplementedException()
    
    def get_col_pos(self, col):
        raise UnimplementedException()
    
    def get_row_count(self):
        raise UnimplementedException()
        
    def get_col_count(self):
        raise UnimplementedException()
    
    def get_row(self, row_name):
        raise UnimplementedException()

    def get_col(self, col_name):
        raise UnimplementedException()


def cg_new(type_str):
    """
    cg_new takes a type string and creates a new object from the 
    class named, it uses an internally defined map to find all
    official CGData data types. So if a 'genomicMatrix' is requested
    a CGData.GenomicMatrix.GenomicMatrix is initialized.
    
    type_str -- A string name of a CGData type, ie 'genomicMatrix'
    """
    mod_name, cls_name = OBJECT_MAP[type_str]
    module = __import__(mod_name, globals(), locals(), [ cls_name ])
    cls = getattr(module, cls_name)
    out = cls()
    return out

def load(path, zip=None):
    """
    load is a the automatic CGData loading function. There has to 
    be a '.json' file for this function to work. It inspects the 
    '.json' file and uses the 'type' field to determine the 
    appropriate object loader to use. The object is created 
    (using the cg_new function) and the 'read' method is passed
    a handle to the data file. If the 'zip' parameter is not None, 
    then it is used as the path to a zipfile, and the path parameter 
    is used as an path inside the zip file to the object data
    
    path -- path to file (in file system space if zip is None, otherwise
    it is the location in the zip file)
    zip -- path to zip file (None by default)
    """
    if not path.endswith(".json"):
        path = path + ".json"

    data_path = re.sub(r'.json$', '', path)

    try:
        handle = open(path, 'rU')
        meta = json.loads(handle.read())
    except IOError:
        raise FormatException("Meta-info (%s) file not found" % (path))

    # Throw away empty values
    meta = dict((k, v) for k, v in meta.iteritems() if v != None)

    if meta['cgdata']['type'] in OBJECT_MAP:
        out = cg_new(meta['cgdata']['type'])
        out.update( meta )
        out.path = data_path
        out.load(data_path)
        return out
    else:
        raise FormatException("%s class not found" % (meta['cgdata']['type']))


def light_load(path, zip=None):
    if not path.endswith(".json"):
        path = path + ".json"

    data_path = re.sub(r'.json$', '', path)

    if zip is None:
        try:
            handle = open(path, 'rU')
            meta = json.loads(handle.read())
        except IOError:
            raise FormatException("Meta-info (%s) file not found" % (path))
    else:
        z = ZipFile(zip)
        handle = z.open(path,'rU')
        meta = json.loads(handle.read())
        handle.close()
        z.close()

    # Throw away empty values
    meta = dict((k, v) for k, v in meta.iteritems() if v != None)

    if meta['cgdata']['type'] in OBJECT_MAP:
        out = cg_new(meta['cgdata']['type'])
        out.update( meta )
        out.path = data_path
        out.zip = zip
        out.light_mode = True
        return out
    else:
        raise FormatException("%s class not found" % (meta['cgdata']['type']))

global LOG_LEVEL
LOG_LEVEL = 2

def info(eStr):
    if LOG_LEVEL < 2:
        sys.stderr.write("LOG: %s\n" % (eStr))
        #errorLogHandle.write("LOG: %s\n" % (eStr))

def debug(eStr):
    if LOG_LEVEL < 1:
        sys.stderr.write("DEBUG: %s\n" % (eStr))
        #errorLogHandle.write("LOG: %s\n" % (eStr))

def warn(eStr):
    if LOG_LEVEL < 3:
        sys.stderr.write("WARNING: %s\n" % (eStr))
        #errorLogHandle.write("WARNING: %s\n" % (eStr))


def error(eStr):
    sys.stderr.write("ERROR: %s\n" % (eStr))
    #errorLogHandle.write("ERROR: %s\n" % (eStr))

