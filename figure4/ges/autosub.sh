#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for const in 0.125 0.25 0.5 1 2 4 8 16 32 64 70 80 90 100
	do
		file=`basename $file`
		echo "R --no-save --args ../data/"${file} ${const}" < ges.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.${const}.sh
		rm tmp.sh
		sbatch scripts/${file}.${const}.sh
	done
done
