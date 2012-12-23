
import CGData
import CGData.BaseTable

class GenomicSegment(CGData.BaseTable.BaseTable):

    __format__ = {
            "name" : "genomicSegment",
            "type" : "type",
            "form" : "table",
            "columnOrder" : [
                "id",
                "chrom",
                "chrom_start",
                "chrom_end",
                "strand",
                "value"
            ],
            "groupKey" : "id",
            "columnDef" : {
                "chrom_start" : { "type" : "int" },
                "chrom_end" :   { "type" : "int" },
                "value" : { "type" : "float" }
            },
            "links" : {
                "assembly" : {},
                "dataSubType" : {}
            } 
        }

    def __init__(self):
        CGData.BaseTable.BaseTable.__init__(self)
