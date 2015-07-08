#change information in icgcLib.py
#change cohortMeta.tsv 
./download.py
./ids.py
./snpEff.py
./snpEffParse.py
./xenaMatrix.py
./rewriteMetadata.py


python /inside/home/jzhu/scripts/mergeGenomicMatrixFiles.py /data/TCGA/icgcFiles/xena/mutGene.icgc /data/TCGA/icgcFiles/xena/mutGene.*.tsv

python /inside/home/jzhu/scripts/mergeMultipleXenaMutation.py /data/TCGA/icgcFiles/xena/simple_somatic_mutation.open.icgc cohort.icgc errorlog /data/TCGA/icgcFiles/xena/simple_somatic_mutation.open.*tsv

python /inside/home/jzhu/scripts/mergeGenomicMatrixFiles.py /data/TCGA/icgcFiles/xena/exp_seq.icgc /data/TCGA/icgcFiles/xena/exp_seq.*.tsv

python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/clinical.icgc /data/TCGA/icgcFiles/xena/clinical*tsv

python cohort.py /data/TCGA/icgcFiles/xena/phenotype.icgc /data/TCGA/icgcFiles/xena/clinical*.tsv