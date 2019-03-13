if [ "$#" -ne 1 ]; then
    echo "usage: ./HCA_h5_to_xena_small.bash HCAdata_dir"
    echo
    exit
fi

pythonDir="~/xenaH5"

echo $1
inputdir=$1
for dir in $inputdir/*; do
    if [[ -d $dir ]]; then
	input_h5=$dir/filtered_gene_bc_matrices_h5.h5
	output=$dir/filtered_gene_bc_matrices.tsv
	prefix=`basename $dir`
	python $pythonDir/h5_xena.py $input_h5 GRCh38 $output $prefix &
    fi
done
