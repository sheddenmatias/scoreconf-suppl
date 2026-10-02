#!/bin/bash
#SBATCH --job-name=reconf_wmerge   #Job name
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/reconf/Output/merge_width.%j.out

# Pools the width blocks into width_dat.Rds, which the fetch downloads.
# Submitted with --dependency=afterany, so it reports short cells rather than
# refusing to produce anything when a task failed.

date; hostname; pwd

module load R

export RECONF_DIR=/blue/k.ekvall/k.ekvall/Simulations/reconf
export RECONF_WIDTHDIR=$RECONF_DIR/Results/width_blocks
export RECONF_RESDIR=$RECONF_DIR/Results
export RECONF_BLOCK=250

Rscript "$RECONF_DIR/merge_width.R"

date
