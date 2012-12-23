#!/usr/bin/env python

import sys
import CGData.GenomicSegment
import CGData.GeneMap
import CGData.RefGene
import CGData.GenomicMatrix

if __name__ == "__main__":

    sg = CGData.GenomicSegment.GenomicSegment()
    sg.load( sys.argv[1]  )

    rg = CGData.RefGene.RefGene()
    rg.load( sys.argv[2] )

    pm = CGData.GeneMap.ProbeMapper('b')

    out = CGData.GeneMap.genomicSegment2MatrixNorm(sg,rg,pm)
    out.write(sys.stdout)

