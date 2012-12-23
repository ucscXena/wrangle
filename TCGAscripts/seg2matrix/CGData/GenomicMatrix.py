
import csv
import CGData
import CGData.BaseMatrix

class GenomicMatrix(CGData.BaseMatrix.BaseMatrix):

    __format__ = {
            "name" : "genomicMatrix",
            "type" : "type",
            "form" : "matrix",
            "rowType" : "probeMap",
            "colType" : "idMap",
            "valueType" : "float",
            "nullString" : "NA",
            "links" : {
                "dataSubType" : {}
            }
        }    

    def __init__(self):
        CGData.BaseMatrix.BaseMatrix.__init__(self)

    def init_blank(self, cols, rows):
        super(GenomicMatrix, self).init_blank(cols=cols,rows=rows)
        if 'cgdata' not in self:
            self['cgdata'] = { 'type' : 'genomicMatrix' }

    def get_probe_list(self):
        return self.get_row_list()

    def get_sample_list(self):
        return self.get_col_list()
    
    def get_data_subtype(self):
        return self.get('cgdata', {}).get('dataSubType', None)

