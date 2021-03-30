#!/usr/bin/env bash

for file in `ls data/*_neig_?.rda`
do
	for lambda in 0.001 0.01 0.1
	do
		file=`basename $file`
		echo "R --no-save --args data/"${file} ${lambda}" < methods.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.${lambda}.methods.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.methods.sh
	done
done
