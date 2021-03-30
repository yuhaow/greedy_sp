#!/usr/bin/env bash

for file in `ls data/*_neig_?.rda`
do
	for lambda in 0.001 0.01 0.1
	do
		file=`basename $file`
		echo "R --no-save --args data/"${file} ${lambda}" < sp_d5.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.${lambda}.sp.sh
		rm tmp.sh
		sbatch scripts/${file}.${lambda}.sp.sh
	done
done
