if [ "$#" -ne 1 ]; then
    echo "usage: ./ensemble.gene.transcript.probemap.bash ensemble_downloaded_gtf"
    echo
    exit
fi

echo $1

file=$1
root=$(basename $file)

genefile=$root.gene
transcriptfile=$root.transcript

geneProbemap=$genefile.probemap
transcriptProbemap=$transcriptfile.probemap

awk '$3 == "gene"' $file > $genefile
awk '$3 == "transcript"' $file > $transcriptfile

printf '%b\t%b\t%b\t%b\t%b\t%b\n' "id" "gene" "chrom" "chromStart" "chromEnd" "strand" > $geneProbemap
cut -f 1,4,5,7,9 $genefile | \
awk 'BEGIN{FS="\t"; OFS="\t"} 
	{split($5,a,";");split(a[1],b," ");split(a[3],c," "); split(a[2],d," "); 
	print substr(b[2],2,length(b[2])-2) "." substr(d[2],2,length(d[2])-2), substr(c[2],2, length(c[2])-2),$1,$2,$3,$4}' >> $geneProbemap

printf '%b\t%b\t%b\t%b\t%b\t%b\n' "id" "gene" "chrom" "chromStart" "chromEnd" "strand" > $transcriptProbemap
cut -f 1,4,5,7,9 $transcriptfile | \
awk 'BEGIN{FS="\t"; OFS="\t"} 
	{split($5,a,";");split(a[3],b," ");split(a[5],c," "); split(a[4],d," "); 
	print substr(b[2],2,length(b[2])-2) "." substr(d[2],2,length(d[2])-2), substr(c[2],2, length(c[2])-2),$1,$2,$3,$4}' >> $transcriptProbemap
