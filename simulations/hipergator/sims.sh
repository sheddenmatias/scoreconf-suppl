#!/bin/bash
#SBATCH --job-name=reconf_cov   #Job name
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --array=1-20
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/reconf/Output/r_job.%j_%a.out

# Coverage sweeps for Figure 1: 20 cells across two sweeps, one
# cell per array task. simulations/code/cells.R holds the cell table and
# run_cell.R does the work of one task.
#
#   sbatch sims.sh
#
# is the whole submission -- the pooling job queues itself, below.
#
# Resubmitting the same array is how a run is finished, not a repair: blocks
# already written are skipped, so only what timed out or failed is redone.
#
# --array=1-20   one task per cell, which is the finest granularity available
#                and so the best defence against the spread in cell cost. Add
#                %10 to cap concurrency if the burst allocation cannot hold 20
#                tasks at once; the array still completes, in more passes.
#                Keep equal to the number of cells: run_cell.R strides, so a
#                smaller array puts two cells on the low-numbered tasks.
# --time         twelve hours, raised from four when the replication count
#                went from 2000 to 10000. Short jobs backfill into idle time on
#                a low-priority burst QOS where long ones wait for a clean
#                window, so keep this as short as the slowest cell allows; the
#                block saving in run_cell.R makes being cut off cheap.

#Record the time and compute node the job ran on
date; hostname; pwd

#Use modules to load the environment for R
module load R

export RECONF_DIR=/blue/k.ekvall/k.ekvall/Simulations/reconf
export RECONF_CELLDIR=$RECONF_DIR/Results/cells   # one file per block
export RECONF_BLOCK=250

mkdir -p "$RECONF_CELLDIR"

# One task queues the pooling job, held until the whole array finishes, so that
# submitting this script is the only command needed. It goes before the work
# below, not after, so the pooling job is queued even if this task then fails --
# its report of which cells are short is exactly what a failed run needs.
# Resubmitting to fill gaps queues a fresh pooling job the same way.
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    sbatch --dependency=afterany:$SLURM_ARRAY_JOB_ID \
           "$RECONF_DIR/merge.sh"
fi

#Run R script
Rscript "$RECONF_DIR/run_cell.R"

date
