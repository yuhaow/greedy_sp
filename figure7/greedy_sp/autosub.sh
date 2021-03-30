#!/bin/sh

for file in `ls ./greedy_sp_data/*.rda`
do

	file=`basename $file`
	echo "R --no-save --args ./greedy_sp_data/"${file}" < sp.R" > tmp.sh
	cat template.sh tmp.sh > scripts/${file}.sh
	rm tmp.sh
	sbatch scripts/${file}.sh

done
