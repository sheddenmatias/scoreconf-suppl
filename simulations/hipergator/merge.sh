#!/bin/bash
#SBATCH --job-name=reconf_merge   #Job name
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/reconf/Output/merge.%j.out

# Pools the per-block files into cov_dat.Rds and covrho_dat.Rds -- the files
# fig_routes.R reads.
#
# Submit with --dependency=afterany rather than afterok, so it runs even when
# some tasks failed: it names the short cells instead of refusing to produce
# anything, which is what tells you which array indices to resubmit.

date; hostname; pwd

module load R

export RECONF_DIR=/blue/k.ekvall/k.ekvall/Simulations/reconf
export RECONF_CELLDIR=$RECONF_DIR/Results/cells
export RECONF_RESDIR=$RECONF_DIR/Results
export RECONF_BLOCK=250

Rscript "$RECONF_DIR/merge_cells.R"

date
