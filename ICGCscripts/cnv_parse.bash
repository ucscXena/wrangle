#!/bin/bash

cancer='*'
#cancer='MELA-AU'

for file in copy_number_somatic_mutation.$cancer.tsv ; do
    output=${file/.tsv/.HiSeq.sp.xena}
    echo $output
    grep "Illumina HiSeq" $file | cut -f 3,10,12,13,14 | sed  '/\t\t/d' | awk  'BEGIN{OFS="\t"} {print $1,$3,$4,$5,$2}' | sed -e 's/[-]\?[0-9]\+.[0-9]\+E-[0-9]\+\b/0/' > $output &

    output=${file/.tsv/.SNP6.sp.xena}
    echo $output
    grep "Affymetrix Genome-Wide Human SNP Array 6.0" $file | cut -f 3,10,12,13,14 | sed  '/\t\t/d' | awk  'BEGIN{OFS="\t"} {print $1,$3,$4,$5,$2}' | sed -e 's/[-]\?[0-9]\+.[0-9]\+E-[0-9]\+\b/0/' > $output &
done
