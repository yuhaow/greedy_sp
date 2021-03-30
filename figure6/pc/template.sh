#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 12:00:00     # 12 hours
#SBATCH -p newnodes     # partition name
#SBATCH -J high-dim-PC  # sensible name for the job

. /etc/profile.d/modules.sh
module add /cm/shared/modulefiles/engaging/gmp/6.1.1

