
import CGData.BaseTable


class ProbeMap(CGData.BaseTable.BaseTable):
    __format__ =  {
            "name" : "probeMap",
            "form" : "table",
            "columnOrder" : [
                "probe",
                "aliases",
                "chrom",
                "chrom_start",
                "chrom_end",
                "strand"
            ],
            "primaryKey" : "probe",
            "columnDef" : {
				"chrom_start" : { "type" : "int", "index" : 1 },
				"chrom_end" : { "type" : "int", "index" : 1 }
            }
        }
        
    def __init__(self):
        CGData.BaseTable.BaseTable.__init__(self)

