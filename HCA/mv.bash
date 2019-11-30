#!/bin/bash

root=/mnt/efsSingleCell/data/HCA
xenaroot=/mnt/efsSingleCell/xena/files/HCA

cd $root

#for dir in $(ls -d */); do
for dir in KidneySingleCellAtlas; do
    echo $dir
    ls $root/$dir/matrix/*.tsv $root/$dir/matrix/*.tsv.json
    ls -l $xenaroot/$dir
    mv $root/$dir/matrix/*.tsv $xenaroot/$dir/
    mv $root/$dir/matrix/*.tsv.json $xenaroot/$dir/
done
