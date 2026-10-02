#!/bin/bash
#SBATCH --job-name=reconf_test   #Job name
#SBATCH --mail-type=FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --array=1-3
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/reconf/Output/test.%j_%a.out

# Submit this before the full array. Three tasks take every third cell, so
# between them they touch all 20 coverage cells, at ten replications
# each. A few minutes, into a throwaway directory, and it cannot touch the real
# results. It exercises the packages, the paths, the array index and the block
# writing. If the three logs end with "task N done", sims.sh will run.

date; hostname; pwd

module load R

export RECONF_DIR=/blue/k.ekvall/k.ekvall/Simulations/reconf
export RECONF_CELLDIR=$RECONF_DIR/Results/test_cells
export RECONF_NSIM=10
export RECONF_BLOCK=5

mkdir -p "$RECONF_CELLDIR"

Rscript "$RECONF_DIR/run_cell.R"

date
