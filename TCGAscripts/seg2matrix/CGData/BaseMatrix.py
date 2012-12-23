
import csv
import CGData
import math
from copy import copy
try:
    import numpy
except ImportError:
    numpy = None

class BaseMatrix(CGData.CGDataMatrixObject):
    """
    Core matrix class. Implements data matrix using numpy or native python objects
    depending up avaliblity and user request
    """
    corner_name = "#"
    element_type = str
    null_type = None
    def __init__(self,type=str):
        CGData.CGDataMatrixObject.__init__(self)
        self.free()
        if 'cgformat' in self and 'valueType' in self['cgformat']:
            if self['cgformat']["valueType"] == 'float':
                self.element_type = float
        else:
            self.element_type = type

    def free(self):
        self.col_map = {}
        self.row_map = {}    
        self.matrix = None
    
    def init_blank(self, cols, rows, skip_numpy=False):
        """
        Initlize matrix with NA (or nan) values using row/column names
        provided by user. User can also force usage of native python objects
        (which is useful for string based matrices, and numpy matrices fix cel string length)
        """
        if numpy is not None and not skip_numpy:
            self.matrix = numpy.matrix( numpy.zeros( (len(rows), len(cols)),  dtype=self.element_type) )
            self.matrix.fill( numpy.nan )
        else:
            self.matrix = []
            for i in range(len(rows)):
				self.matrix.append([self.null_type]*len(cols))
        for i, c in enumerate(cols):
            self.col_map[c] = i
        for i, r in enumerate(rows):
            self.row_map[r] = i
        self.loaded = True

    def read(self, handle, skip_vals=False):
        self.col_map = {}
        self.row_map = {}    
        pos_hash = None

        if numpy is not None:
            #txtMatrix = numpy.loadtxt(handle, delimiter="\t", comments="%%%%%%%%%%%%%%", dtype=str)
            t = []
            for line in handle:
                t.append(line.replace("\n", "").split("\t"))
            txtMatrix = numpy.array(t)
            del t
            if self.element_type == float:
                txtMatrix[ txtMatrix=="NA" ] = 'nan'
                txtMatrix[ txtMatrix=="null" ] = 'nan'
                self.matrix = numpy.matrix( numpy.zeros( (txtMatrix.shape[0]-1, txtMatrix.shape[1]-1) ) )
                self.matrix.fill(numpy.nan)
                for i in range(self.matrix.shape[0]):
                    for j in range(self.matrix.shape[1]):
                        try:
                            self.matrix[i,j] = self.element_type(txtMatrix[i+1,j+1])
                        except ValueError:
                            pass
            else:
                self.matrix = numpy.matrix(txtMatrix[1:,1:], dtype=self.element_type)
                
            for i, col in enumerate( txtMatrix[0,1:] ):
                self.col_map[col] = i
            for i, row in enumerate( txtMatrix[1:,0] ):
                self.row_map[row] = i
        else:
            self.matrix = []
            for row in csv.reader(handle, delimiter="\t"):
                if pos_hash is None:
                    pos_hash = {}
                    pos = 0
                    for name in row[1:]:
                        i = 1
                        orig_name = name
                        while name in pos_hash:
                            name = orig_name + "#" + str(i)
                            i += 1
                        pos_hash[name] = pos
                        pos += 1
                else:
                    newRow = []
                    if not skip_vals:                    
                        newRow = [self.null_type] * (len(pos_hash))
                        for col in pos_hash:
                            i = pos_hash[col] + 1
                            if row[i] != 'NA' and row[i] != 'null' and row[i] != 'NONE' and row[i] != "N/A" and len(row[i]):
                                newRow[i - 1] = self.element_type(row[i])
                    self.row_map[row[0]] = len(self.matrix)
                    self.matrix.append(newRow)

            self.col_map = {}
            for col in pos_hash:
                self.col_map[col] = pos_hash[col]
        self.loaded = True

    def write(self, handle, missing='NA'):
        write = csv.writer(handle, delimiter="\t", lineterminator='\n')
        col_list = self.get_col_list()
        
        write.writerow([self.corner_name] + col_list)
        for rowName in self.row_map:
            out = [rowName]
            row = self.get_row(rowName)
            for col in col_list:
                val = row[self.col_map[col]]
                if val == self.null_type or val is None or (type(val)==float and math.isnan(val)):
                    val = missing
                out.append(val)
            write.writerow(out)
    
    def read_keyset(self, handle, key_predicate):
        if key_predicate == "rowKeySrc":
            reader = csv.reader( handle, delimiter="\t")
            head = None
            for row in reader:
                if head is None:
                    head = row
                else:
                    yield row[0]
        
        if key_predicate=="columnKeySrc":
            reader = csv.reader( handle, delimiter="\t")
            head = None
            for row in reader:
                for col in row[1:]:
                    yield col
                break
                
    def get_col_namespace(self):
        """
        Return the name of the column namespace
        """
        return self.get("colNamespace", None)

    def get_row_namespace(self):
        """
        Return the name of the row namespace
        """
        return self.get("rowNamespace", None)
        
    def get_col_list(self):
        """
        Returns names of columns
        """
        if not self.loaded:
            self.load( )
        out = self.col_map.keys()
        out.sort( lambda x,y: self.col_map[x]-self.col_map[y])
        return out 
        
    def get_row_list(self):
        """
        Returns names of rows
        """
        out = self.row_map.keys()
        out.sort( lambda x,y: self.row_map[x]-self.row_map[y])
        return out 
    
    def get_row_pos(self, row):
        return self.row_map[row]
    
    def get_col_pos(self, col):
        return self.col_map[col]
    
    def get_row_count(self):
        return len(self.row_map)
        
    def get_col_count(self):
        return len(self.col_map)
    
    def get_row_map(self):
        return copy(self.row_map)

    def get_col_map(self):
        return copy(self.col_map)
    
    def get_shape(self):
        return len(self.row_map), len(self.col_map)
    
    def get_row(self, row_name):
        if not self.loaded:
            self.load( )
        if isinstance(self.matrix, list):
            return self.matrix[ self.row_map[row_name] ]
        else:
            return self.matrix[ self.row_map[row_name] ].tolist()[0]

    def get_col(self, col_name):
        if not self.loaded:
            self.load( )
        if isinstance(self.matrix, list):
            out = []
            for row_name in self.get_row_list():
                out.append( self.get_val(col_name, row_name) )
            return out
        else:
            return self.matrix[:,self.col_map[col_name]].reshape(-1).tolist()[0]
    
    def get_val(self, col_name, row_name):
        """
        Get cell value based on row and column names
        """
        if isinstance(self.matrix, list):
            return self.matrix[self.row_map[row_name]][self.col_map[col_name]]
        return self.matrix[self.row_map[row_name],self.col_map[col_name]]
    
    def set_val(self, col_name, row_name, value):
        """
        Set cell value based on row and column names
        """
        if isinstance(self.matrix, list):
            self.matrix[self.row_map[row_name]][self.col_map[col_name]] = value
        else:
            self.matrix[self.row_map[row_name],self.col_map[col_name]] = value

    def write_gct(self, handle, missing=''):
        write = csv.writer(handle, delimiter="\t", lineterminator='\n')
        cols = self.get_col_list()
        write.writerow(["#1.2"])
        write.writerow([len(self.get_row_list()), len(self.get_col_list())])
        write.writerow(["NAME", "Description"] + cols)
        for row in self.get_row_list():
            out = [row, row]
            for col in cols:
                val = self.get_val(row_name=row, col_name=col)
                if val is None:
                    val = missing
                out.append(val)
            write.writerow(out)


