# Simulations on HiPerGator

Three of the result files in `results/` were produced on the University of
Florida's HiPerGator cluster, with SLURM, in August 2026:

- `cov_dat.Rds` and `covrho_dat.Rds`, for Figure 1;
- `width_dat.Rds`, for Figure 2.

This file records what was run and how to run it again. The simulation
for Figure 3 was run locally with `code/boundary_k_sims.R`. 

## Paper Run

| Job | Cells used in the paper | Replications | Produces |
|---|---|---|---|
| `hipergator/sims.sh` | 20 coverage cells: sweeps `cov` and `covrho` in `code/cells.R` | 10,000 per cell, in 40 blocks of 250 | `cov_dat.Rds`, `covrho_dat.Rds` |
| `hipergator/width.sh` | 5 width cells in `code/width_engine.R` | 1,000 per cell, in 4 blocks of 250 | `width_dat.Rds` |

The jobs ran under these settings:

- **Software:** R from the cluster's `R` module (4.6.1 in August 2026), and
  `reconf` installed from GitHub at commit
  `164d7bd377f053aca3177279e71487507e54834d`.
- **Array tasks:** each task ran one cell.
  - Coverage tasks requested 8 cores, 32 GB and 12 hours.
  - Width tasks requested 8 cores, 16 GB and 4 hours.

Note: the saved `cov_dat.Rds` and `width_dat.Rds` contain some extra rows from 
the original run that the figure scripts ignore.

## Running it again on HiPerGator

The scripts are set up for the account `k.ekvall`. To run
under another account, or on another SLURM cluster, change:

| File | Lines to change |
|---|---|
| `sims.sh`, `width.sh`, `merge.sh`, `merge_width.sh`, `test_one.sh` | `--account`, `--qos`, `--mail-user`, the path in `--output`, and `RECONF_DIR` |
| `hpg.sh` | `REMOTE` (user and host) and `BASE` (the cluster path) |
| This file | The user and path in the upload commands below |

With the same code and the same `reconf` commit, a rerun reproduces the
paper's cells. Newer versions of the other packages may change results.

### 1. Upload

The contents of `simulations/code/` and the job scripts in
`simulations/hipergator/` go into one directory. From the repository
root:

```bash
rsync -av simulations/code/ \
          k.ekvall@hpg.rc.ufl.edu:/blue/k.ekvall/k.ekvall/Simulations/reconf/
rsync -av simulations/hipergator/*.sh \
          k.ekvall@hpg.rc.ufl.edu:/blue/k.ekvall/k.ekvall/Simulations/reconf/
ssh k.ekvall@hpg.rc.ufl.edu \
    mkdir -p /blue/k.ekvall/k.ekvall/Simulations/reconf/{Output,Results/cells,Results/width_blocks}
```

### 2. Install the packages

On the cluster, start an interactive session on the development partition
(`srundev --time=01:00:00 --cpus-per-task=4 --mem=16gb`), run
`module load R`, start R, and run:

```r
pkgs <- c("Matrix","lme4","trust","assertthat","generics","Rcpp","RcppEigen",
          "nloptr","alabama","dplyr","here","foreach","doParallel","doRNG","remotes")
install.packages(pkgs, repos = "https://cloud.r-project.org", Ncpus = 4)

remotes::install_github("koekvall/reconf",
                        ref = "164d7bd377f053aca3177279e71487507e54834d")
packageDescription("reconf")$RemoteSha   # must print the same SHA
```

### 3. Clear old blocks, test, then submit

Blocks on disk are matched by cell key only. If blocks from an earlier run
are still in `Results/cells` or `Results/width_blocks`, a new submission
skips those cells and re-pools the old results. Before rerunning from
scratch, delete both folders and recreate them:

```bash
cd /blue/k.ekvall/k.ekvall/Simulations/reconf
rm -rf Results/cells Results/width_blocks && mkdir -p Results/cells Results/width_blocks
```

Then:

```bash
sbatch /blue/k.ekvall/k.ekvall/Simulations/reconf/test_one.sh   # 3 tasks, 10 replications per cell
sbatch /blue/k.ekvall/k.ekvall/Simulations/reconf/sims.sh       # coverage, 20 tasks
sbatch /blue/k.ekvall/k.ekvall/Simulations/reconf/width.sh      # width, 5 tasks
```

- `test_one.sh` writes to `Results/test_cells`, so it cannot touch real
  results. Its three logs (`Output/test.*.out`) should each end with
  "task N done".
- `sims.sh` and `width.sh` each queue their own pooling job, so nothing else
  needs to be submitted.

### 4. Collect

- **Check the pooling logs** (`Output/merge.*.out` and
  `Output/merge_width.*.out`). They name any cell that is short of its 40 or
  4 blocks.
- **If cells are short,** submit the same script again. Only the missing
  blocks run, and a new pooling job is queued.
- **When every cell is complete,** download `cov_dat.Rds`, `covrho_dat.Rds`
  and `width_dat.Rds` from `Results/` into the repository's `results/`.

### 5. Draw the figures

On your own computer, from the repository root:

```bash
Rscript simulations/fig_routes.R    # Figure 1: figures/fig_cov_routes.pdf
Rscript simulations/fig_width.R     # Figure 2: figures/fig_width_dists.pdf
```

The figure scripts need `tidyverse` with `ggplot2` 3.5.0 or later,
`patchwork`, `ggthemes` and `here`.

## Using `hpg.sh`

`hipergator/hpg.sh` handles upload, submission, status check and download
from a laptop. `--fresh` deletes the old blocks before submitting, as in
step 3.

```bash
simulations/hipergator/hpg.sh test
simulations/hipergator/hpg.sh run --fresh
simulations/hipergator/hpg.sh width --fresh
simulations/hipergator/hpg.sh status
simulations/hipergator/hpg.sh fetch    # downloads every .Rds in Results/ into results/
```
