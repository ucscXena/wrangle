
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
    'sampleMap': ('CGData.SampleMap', 'SampleMap'),
    'clinicalMatrix': ('CGData.ClinicalMatrix', 'ClinicalMatrix'),
    'dataSubType': ('CGData.DataSubType', 'DataSubType'),
    'trackGenomic': ('CGData.TrackGenomic', 'TrackGenomic'),
    'trackClinical': ('CGData.TrackClinical', 'TrackClinical'),
    'assembly': ('CGData.Assembly', 'Assembly'),
    'clinicalFeature': ('CGData.ClinicalFeature', 'ClinicalFeature'),
    'refGene' : ('CGData.RefGene', 'RefGene')
}

MERGE_OBJECTS = [ 'trackClinical', 'trackGenomic' ]

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

class CGGroupMember(object):
    pass

class CGGroupBase(object):

    DATA_FORM = None

    def __init__(self, group_name):
        self.members = {}
        self.name = group_name
    
    def __setitem__(self, name, item):
        self.members[ name ] = item
    
    def __getitem__(self, name):
        return self.members[ name ]
    
    def put(self, obj):
        self.members[ obj.get_name() ] = obj
    
    def is_link_ready(self):
        for name in self.members:
            if not self.members[name].is_link_ready():
                return False
        return True
    
    def get_name(self):
        return self.name
    
    def unload(self):
        for name in self.members:
            self.members[name].unload()
    
    def lookup(self, **kw):
        for elem in self.members:
            found = True
            obj = self.members[ elem ]
            for key in kw:
                if obj.get( key, None ) != kw[key] and obj.get( ":" + key, None ) != kw[key]:
                    found = False
            if found:
                return obj
                    
    
    def get_link_map(self):
        out = {}
        for name in self.members:
            lmap = self.members[ name ].get_link_map()
            for ltype in lmap:
                if ltype not in out:
                    out[ ltype ] = []
                for lname in lmap[ltype]:
                    if lname not in out[ltype]:
                        out[ltype].append( lname )
        return out
    
class UnimplementedException(Exception):
    def __init__(self, str):
        Exception.__init__(self, str)

class CGObjectBase(dict):
    """
    This is the base object for CGData loadable objects.
    The methods covered in the base case cover usage meta-information
    loading/unloading and manipulation as well as zip (cgz) file access.
    """
    def __init__(self):
        self.path = None
        self.zip = None
        self.light_mode = False
        super(CGObjectBase,self).__init__()

    # XXX There are no less than three different code paths for
    # loading the json data, that get hit at different points during
    # the compile. This is really messed up. There's this routine, plus
    # the load and light_load methods, below.
    def load(self, path=None, **kw):
        if path is None and self.path is not None:
            path = self.path
        if path is None:
            raise OSError( "Path not defined" ) 
        
        if self.zip is None:
            dhandle = open(path)
            self.read(dhandle, **kw)
            dhandle.close()
        else:
            z = ZipFile(self.zip)
            dhandle = z.open(self.path)
            self.read(dhandle, **kw)
            dhandle.close()
            z.close()
            
        self.path = path
        if (os.path.exists(path + ".json")):
            mhandle = open(path + ".json")
            meta = json.loads(mhandle.read())
            meta = dict((k, v) for k, v in meta.iteritems() if v != None)
            self.update(meta)
            mhandle.close()

    def unload(self):
        pass

    def is_link_ready(self):
        return True
    
    def store(self, path=None):
        if path is None and self.path is not None:
            path = self.path
        if path is None:
            raise OSError( "Path not defined" ) 
        mHandle = open(path + ".json", "w")
        mHandle.write(json.dumps(self))
        mHandle.close()
        if not self.light_mode:
            self.path = path
            dhandle = open(path, "w")
            self.write(dhandle)
            dhandle.close()            
    
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
    
    def is_group_member(self):
        if 'group' in self:
            return True
        return False
    
    def get_group(self):
        return self.get( 'group', self.get('name', None))

    def get_name(self):
        return self.get( 'name', None )
    
    def get_link_map(self):
        out = {}
        for key in self:
            if key.startswith(':'):
                if isinstance( self[ key ], list ):
                    out[ key[1:] ] = self[ key ]
                elif self[ key ] is not None:
                    out[ key[1:] ] = [ self[ key ] ]
        return out

    def add_history(self, desc):
        if not 'history' in self:
            self[ 'history' ] = []
        self[ 'history' ].append( desc )

    

class CGMergeObject(object):
    
    typeSet = {}
    
    def __init__(self):
        self.members = {}
    
    def merge(self, **kw):
        self.members = kw
    
    def __iter__(self):
        return self.members.keys().__iter__()
    
    def __getitem__(self, item):
        return self.members[item]

    def unload(self):
        pass
        
    def sql_pass(self, id_table, method):
        for t in self.members:
            if hasattr(self.members[t], "gen_sql_" + method):
                f = getattr(self.members[t], "gen_sql_" + method)
                for line in f(id_table):
                    yield line



class CGDataSetObject(CGObjectBase):
    
    def __init__(self):
        CGObjectBase.__init__(self)


class CGDataMatrixObject(CGObjectBase):
        
    def __init__(self):
        CGObjectBase.__init__(self)


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
        handle = open(path)
        meta = json.loads(handle.read())
    except IOError:
        raise FormatException("Meta-info (%s) file not found" % (path))

    # Throw away empty values
    meta = dict((k, v) for k, v in meta.iteritems() if v != None)

    if meta['type'] in OBJECT_MAP:
        out = cg_new(meta['type'])
        out.update( meta )
        out.path = data_path
        out.load(data_path)
        return out
    else:
        raise FormatException("%s class not found" % (meta['type']))


def light_load(path, zip=None):
    if not path.endswith(".json"):
        path = path + ".json"

    data_path = re.sub(r'.json$', '', path)

    if zip is None:
        try:
            handle = open(path)
            meta = json.loads(handle.read())
        except IOError:
            raise FormatException("Meta-info (%s) file not found" % (path))
    else:
        z = ZipFile(zip)
        handle = z.open(path)
        meta = json.loads(handle.read())
        handle.close()
        z.close()

    # Throw away empty values
    meta = dict((k, v) for k, v in meta.iteritems() if v != None)

    if meta['type'] in OBJECT_MAP:
        out = cg_new(meta['type'])
        out.update( meta )
        out.path = data_path
        out.zip = zip
        out.light_mode = True
        return out
    else:
        raise FormatException("%s class not found" % (meta['type']))

global LOG_LEVEL
LOG_LEVEL = 2

def log(eStr):
    if LOG_LEVEL < 2:
        sys.stderr.write("LOG: %s\n" % (eStr))
        #errorLogHandle.write("LOG: %s\n" % (eStr))


def warn(eStr):
    if LOG_LEVEL < 1:
        sys.stderr.write("WARNING: %s\n" % (eStr))
        #errorLogHandle.write("WARNING: %s\n" % (eStr))


def error(eStr):
    sys.stderr.write("ERROR: %s\n" % (eStr))
    #errorLogHandle.write("ERROR: %s\n" % (eStr))


#####################


TABLE = "table"
MATRIX = "matrix"

class Column(object):
    def __init__(self, name, type, primary_key=False):
        self.name = name
        self.type = type
        self.primary_key = primary_key

