#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	#for alpha in 0.1 0.2
	for alpha in 0.25 0.3 0.4 0.5 0.6
	do
		file=`basename $file`
		echo "R --no-save --args ../data/"${file} ${alpha}" < sp.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.${alpha}.sh
		rm tmp.sh
		sbatch scripts/${file}.${alpha}.sh
	done
done
