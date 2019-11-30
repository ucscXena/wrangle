#!/bin/bash

root=/mnt/efsSingleCell/data/HCA

#for dir in SingleCellLiverLandscape KidneySingleCellAtlas; do
for dir in KidneySingleCellAtlas; do
    datadir=$root/$dir/matrix
    echo $datadir
    ls $datadir
    python ~/wrangle/HCA/csv_run.py $datadir 1 1
done
