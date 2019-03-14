if [ "$#" -ne 1 ]; then
    echo "usage: ./HCA_h5_to_xena_small.bash HCAdata_dir"
    echo
    exit
fi

xenaH5Dir=$HOME/xenaH5
HCAcodeDir=$HOME/wrangle/HCA

numRUN=6  # number of simultaneous runs, so not overwhelming the server

inputdir=$1
echo $inputdir
numdir=`ls -1 $inputdir | wc -l`
echo total dir: $numdir
numPerRUN=$(($numdir/$numRUN))
echo total RUN: $numRUN
echo perRUN: $numPerRUN

dirs=($(echo $inputdir/* | tr ' ' '\n'))

for (( i=0; i<$numdir+$numPerRUN-1; i=i+$numPerRUN )) do
  command=''
  for (( j=0; j<$numPerRUN; j++ )) do
     index=$(($i+$j))
     if [[ index -lt numdir ]]; then
	 dir=${dirs[index]}
	 if [[ -d $dir ]]; then
	     echo $index, $dir
	     input_h5=$dir/filtered_gene_bc_matrices_h5.h5
	     output=$dir/filtered_gene_bc_matrices.tsv
	     prefix=`python $HCAcodeDir/cell_suspension_id.py $dir/cell_suspension_0.json`
	     command="$command nice python $xenaH5Dir/h5_xena.py $input_h5 GRCh38 $output $prefix ; "
	 fi
     fi
  done
  if [[ $command != '' ]]; then
      echo $command
      eval $command &
  fi
done
