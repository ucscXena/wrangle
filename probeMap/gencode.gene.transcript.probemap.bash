#!/bin/bash
echo $1

file=$1
root=$(basename $file .gtf)

genefile=$root.gene
transcriptfile=$root.transcript

geneProbemap=$root.gene.probemap
transcriptProbemap=$root.transcript.probemap

awk '$3 == "gene"' $file > $genefile
awk '$3 == "transcript"' $file > $transcriptfile

echo -e 'id\tgene\tchrom\tchromStart\tchromEnd\tstrand' > $geneProbemap
cut -f 1,4,5,7,9 $genefile |awk 'BEGIN{FS="\t"; OFS="\t"}{split($5,a,";");split(a[1],b," ");split(a[3],c," ");print substr(b[2],2,length(b[2])-2),substr(c[2],2, length(c[2])-2),$1,$2,$3,$4}' >> $geneProbemap

echo -e 'id\tgene\tchrom\tchromStart\tchromEnd\tstrand' > $transcriptProbemap
cut -f 1,4,5,7,9 $transcriptfile |awk 'BEGIN{FS="\t"; OFS="\t"}{split($5,a,";");split(a[2],b," ");split(a[4],c," ");print substr(b[2],2,length(b[2])-2),substr(c[2],2, length(c[2])-2),$1,$2,$3,$4}' >> $transcriptProbemap

