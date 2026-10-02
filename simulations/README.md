# Simulations

This folder has the code for the three simulation figures in Section 4.

| Paper figure | Result | Figure script | Simulation |
|---|---|---|---|
| Figure 1 | Coverage for the intercept–slope covariance as (a) the variances go to zero and (b) the correlation goes to one | `fig_routes.R` | coverage sweeps (`code/run_cell.R`) |
| Figure 2 | Widths of the likelihood-ratio and proposed intervals | `fig_width.R` | width comparison (`code/run_width.R`) |
| Figure 3 | Coverage as the number of variance parameters grows | `fig_boundary_k.R` | `code/boundary_k_sims.R` |

## Redraw the paper's figures

The results behind every figure are included in `results/`. From the
repository root:

```bash
Rscript simulations/fig_routes.R        # figures/fig_cov_routes.pdf
Rscript simulations/fig_width.R         # figures/fig_width_dists.pdf
Rscript simulations/fig_boundary_k.R    # figures/fig_cov_boundary_k_range.pdf
```

Each script also prints the numbers quoted in Section 4.

## Rerun a figure's simulations locally

`run_local.R` runs the simulations behind one figure, with the paper's
settings and seeds, and then draws that figure from the new results.

```r
FIGURE <- "routes"   # "routes", "width", "boundary_k" or "all"
NSIM   <- 200        # replications per setting. NA for the paper's
CORES  <- NA         # parallel cores. NA for all cores but one
FRESH  <- FALSE      # TRUE to delete this run's earlier output first
```

Output goes to `local_runs/<figure>_nsim<N>/` in the repository root, in the
subfolders `results/` and `figures/`, so the paper's `results/` and `figures/`
are unchanged.

At the paper's replication counts, `boundary_k` takes about one core-hour.
`routes` and `width` take far longer and were run on a cluster. Start with a
small `NSIM`.

## Run on a SLURM cluster

The coverage sweeps and the width comparison were run on the University of
Florida's HiPerGator cluster with the job scripts in `hipergator/`.
`hipergator/hipergator_info.md` records that run and explains how to run it
again there or on another SLURM cluster.

## Contents

```
simulations/
  README.md             this file
  run_local.R           run one figure's simulations locally and draw it
  fig_routes.R          Figure 1
  fig_width.R           Figure 2
  fig_boundary_k.R      Figure 3
  code/                 simulation code. Runs on a laptop or a cluster
                        (boundary_k_sims.R is the simulation for Figure 3)
  hipergator/           SLURM job scripts, hpg.sh (run HiPerGator sims),
                        and hipergator_info.md (additional info)
```
