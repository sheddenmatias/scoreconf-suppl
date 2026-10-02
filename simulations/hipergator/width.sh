#!/bin/bash
#SBATCH --job-name=reconf_width   #Job name
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --array=1-5
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/reconf/Output/r_width.%j_%a.out

# Width comparison: 5 cells, one per array task; run_width.R does the work
# of one task and width_engine.R holds the cell table. Blocks already written
# are skipped, so resubmitting the same array finishes an interrupted run.

date; hostname; pwd

module load R

export RECONF_DIR=/blue/k.ekvall/k.ekvall/Simulations/reconf
export RECONF_WIDTHDIR=$RECONF_DIR/Results/width_blocks
export RECONF_BLOCK=250

mkdir -p "$RECONF_WIDTHDIR"

# One task queues the pooling job, held until the whole array finishes; see
# sims.sh for why it is queued first and with afterany.
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    sbatch --dependency=afterany:$SLURM_ARRAY_JOB_ID \
           "$RECONF_DIR/merge_width.sh"
fi

Rscript "$RECONF_DIR/run_width.R"

date
