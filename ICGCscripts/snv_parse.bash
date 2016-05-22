#!/bin/bash

for file in simple_somatic_mutation.open.*tsv ; do
    output=${file/.tsv/.xena}
    tmpfile=${file/.tsv/.tmp}
    echo $output
    sed 1d $file |cut -f 4,9,10,11,15,17,26,27,29,30  | grep -v intron_variant |grep -v downstream_gene_variant |grep -v intergenic_region |sort |uniq > $tmpfile &&  python ~/cgDataJing/ICGCscripts/icgcSNVcleanup.py $tmpfile $output && rm -f $tmpfile &
done


