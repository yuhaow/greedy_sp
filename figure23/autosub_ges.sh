#!/usr/bin/env bash

for file in `ls data/*_n_*.rda`
do
	file=`basename $file`
	echo "R --no-save --args data/"${file}" < ges.R" > tmp.sh
	cat template.sh tmp.sh > scripts/${file}.ges.sh
	rm tmp.sh
	sbatch scripts/${file}.ges.sh
done
