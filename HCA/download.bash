#!/bin/bash

#works only with new data

root=/mnt/efsSingleCell/data/HCA

function download {
    dir=$1
    url=$2
    mkdir $root/$dir
    wget $url -O $root/$dir/matrix.zip
    unzip $root/$dir/matrix.zip -d $root/$dir
    rm $root/$dir/matrix.zip
    unzipdir=$(ls $root/$dir/)
    mv $root/$dir/$unzipdir $root/$dir/matrix
}

dir=KidneySingleCellAtlas
projectId=abe1a013-af7a-45ed-8c26-f3793c24a1f4

url=https://data.humancellatlas.org/project-assets/project-matrices/$projectId.csv.zip
download $dir $url

: '
dir=SingleCellLiverLandscape
projectId=4d6f6c96-2a83-43d8-8fe1-0f53bffd4674

dir=pancrease_quake
projectId=cddab57b-6868-4be4-806f-395ed9dd635a

dir=interneuronDevelopment_Lein
projectId=2043c65a-1cf8-4828-a656-9e247d4e64f1

dir=Fetal_Maternal_Interface
projectId=f83165c5-e2ea-4d15-a5cf-33f3550bffde

dir=HPSI_human_cerebral_organoids
projectId=005d611a-14d5-4fbf-846e-571a1f874f70

dir=Tissue_stability
projectId=c4077b3c-5c98-4d26-a614-246d12c2e5d7

dir=HumanTissueTcellActivation
projectId=4a95101c-9ffc-4f30-a809-f04518a23803

dir=WongAdultRetina
projectId=8185730f-4113-40d3-9cc3-929271784c2b

dir=Census_of_Immune_Cells
projectId=cc95ff89-2e68-4a08-a234-480eca21ce79

dir=Human_Hematopoietic_Profiling
projectId=091cf39b-01bc-42e5-9437-f419a66c8a45

dir=KidneySingleCellAtlas
projectId=abe1a013-af7a-45ed-8c26-f3793c24a1f4

dir=HumanColonicMesenchymeIBD
projectId=f8aa201c-4ff1-45a4-890e-840d63459ca2
'
