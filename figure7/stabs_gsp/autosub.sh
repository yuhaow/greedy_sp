#!/bin/sh

echo "R --no-save --args ./greedy_sp_data/Cebpb__excluded_data.rda < stability_selection.R" > tmp.sh
cat template.sh tmp.sh > scripts/${file}.sh
rm tmp.sh
sbatch scripts/${file}.sh
