#change information in icgcLib.py
#change cohortMeta.tsv 
./download.py
./ids.py
./snpEff.py
./snpEffParse.py
./xenaMatrix.py
./rewriteMetadata.py


# copy over all the icgc pan-cancer .json file back and change version
python /inside/home/jzhu/scripts/mergeGenomicMatrixFiles.py /data/TCGA/icgcFiles/xena/mutGene.icgc /data/TCGA/icgcFiles/xena/mutGene.*.tsv

python /inside/home/jzhu/scripts/mergeMultipleXenaMutation.py /data/TCGA/icgcFiles/xena/simple_somatic_mutation.open.icgc cohort.icgc errorlog /data/TCGA/icgcFiles/xena/simple_somatic_mutation.open.*tsv

python /inside/home/jzhu/scripts/mergeGenomicMatrixFiles.py /data/TCGA/icgcFiles/xena/exp_seq.icgc /data/TCGA/icgcFiles/xena/exp_seq.*.tsv

python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/donor.icgc /data/TCGA/icgcFiles/xena/donor.*tsv
python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/donor_exposure.icgc /data/TCGA/icgcFiles/xena/donor_exposure.*tsv
python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/donor_therapy.icgc /data/TCGA/icgcFiles/xena/donor_therapy.*tsv
python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/donor_family.icgc /data/TCGA/icgcFiles/xena/donor_family.*tsv
python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/specimen.icgc /data/TCGA/icgcFiles/xena/specimen.*tsv
python /inside/home/jzhu/scripts/mergeClinicalMatrixFiles.py /data/TCGA/icgcFiles/xena/survival.icgc /data/TCGA/icgcFiles/xena/survival.*tsv

python cohort.py /data/TCGA/icgcFiles/xena/study.icgc /data/TCGA/icgcFiles/xena/donor.*.tsv
