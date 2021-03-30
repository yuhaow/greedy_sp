#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -t 12:00:00     # 12 hours
#SBATCH -p newnodes     # partition name
#SBATCH -J node_10_sp  # sensible name for the job

. /etc/profile.d/modules.sh
module add /cm/shared/modulefiles/engaging/gmp/6.1.1
R=/home/yuhaow/software/R/bin/R

$R --no-save --args ./greedy_sp_data/Cebpb__excluded_data.rda < stability_selection.R
