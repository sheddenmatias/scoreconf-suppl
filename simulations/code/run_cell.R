# Simulate the coverage settings (cells) behind Figure 1. 
# Runs the 20 cells listed in cells.R
# merge_cells.R then pools the output into results/cov_dat.Rds (panel a) and 
# results/covrho_dat.Rds (panel b), which fig_routes.R plots.
#
# On HiPerGator, sims.sh starts one SLURM array task per cell, and each task
# runs the cell whose index matches its array index. Without SLURM, as when
# run_local.R calls it, the script runs every cell in turn.
#
# Each cell's replications are run in blocks of 250, and each block is written
# as soon as it finishes. A task that hits its wall time loses at most one
# block, and resubmitting picks up where it stopped: blocks already on disk are
# skipped. Block b of a cell is seeded with the cell's seed plus b, so a
# resubmitted block reproduces exactly and the pooled result does not depend on
# interruptions.
#
# Environment:
#   SLURM_ARRAY_TASK_ID, SLURM_ARRAY_TASK_COUNT  set by sbatch
#   SLURM_CPUS_PER_TASK, RECONF_CORES            workers for the inner loop
#   RECONF_DIR                                   where the scripts are
#   RECONF_BLOCK                                 replications per block (250)
#   RECONF_SWEEP                                 run only cov, covrho, or both
#   RECONF_NSIM                                  replications per cell (10000)
#   RECONF_CELLDIR                               where the blocks are written
# Leave RECONF_NSIM and RECONF_BLOCK at their defaults to reproduce the paper:
# the seeds are tied to the blocks.

root <- normalizePath(Sys.getenv("RECONF_DIR", unset = getwd()))
setwd(root)

# The scripts sit in the repository's simulations/code/ on a laptop, and
# flat in one directory on a cluster, where nothing but the simulation is being
# run. Helpers are found by name in either layout.
locate <- function(name) {
  for (p in c(file.path(root, name),
              file.path(root, "R", name),
              file.path(root, "simulations", "code", name)))
    if (file.exists(p)) return(p)
  stop("cannot find ", name, " under ", root)
}

# coverage_engine.R reaches its neighbours through here(), which normally finds
# a project root by looking for a .git or .Rproj marker. A batch job starts in
# the home directory and an upload need not carry those markers, so here() is
# defined in the global environment first, where it takes precedence over the
# package of the same name, and resolves by name instead.
here <- function(...) locate(basename(file.path(...)))

source(locate("coverage_engine.R"))
source(locate("cells.R"))

env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(v) || v < 1L) default else v
}

array_id  <- env_int("SLURM_ARRAY_TASK_ID", 1L)
n_arrays  <- env_int("SLURM_ARRAY_TASK_COUNT", 1L)
n_cores   <- env_int("SLURM_CPUS_PER_TASK",
                     env_int("RECONF_CORES", max(1L, parallel::detectCores() - 1L)))
block_sz  <- env_int("RECONF_BLOCK", 250L)

cell_dir <- Sys.getenv("RECONF_CELLDIR",
                       unset = file.path(root, "results", "cells"))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

# RECONF_SWEEP restricts the run to one or more sweeps, comma separated: cov,
# covrho. Unset means all of them, which is what the array wants.
sweeps <- Sys.getenv("RECONF_SWEEP")
idx <- seq_along(CELLS)
if (nzchar(sweeps)) {
  want <- trimws(strsplit(sweeps, ",")[[1]])
  idx <- idx[vapply(CELLS, `[[`, character(1), "sweep") %in% want]
  if (!length(idx)) stop("no cells in sweep(s): ", sweeps)
}

# Task k takes every n_arrays-th cell starting at k. Striding rather than
# taking a contiguous block, so that neighbouring cells of similar cost are
# spread across tasks instead of landing on one.
pos  <- seq_along(idx)
mine <- idx[pos[(pos - array_id) %% n_arrays == 0]]

cat(sprintf("task %d of %d: %d cells, %d cores, blocks of %d\n",
            array_id, n_arrays, length(mine), n_cores, block_sz))

# One cluster for the whole task, reused across cells and blocks: the workers
# load lme4 and reconf on first use and setting one up per block would repeat
# that. Single-threaded BLAS on the workers, since the outer loop already uses
# every core.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
.cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(.cl)
invisible(parallel::clusterEvalQ(.cl, Sys.setenv(VECLIB_MAXIMUM_THREADS = "1",
                                                 OMP_NUM_THREADS = "1")))
parallel::clusterExport(.cl, ls(globalenv()), envir = globalenv())
options(reconf.cluster = .cl)
on.exit({parallel::stopCluster(.cl); options(reconf.cluster = NULL)})

for (ci in mine) {
  cl_spec  <- CELLS[[ci]]
  n_sim    <- env_int("RECONF_NSIM", cl_spec$n_sim)
  n_blocks <- ceiling(n_sim / block_sz)

  for (b in seq_len(n_blocks)) {
    out_file <- file.path(cell_dir, sprintf("%s__b%02d.Rds", cl_spec$key, b))
    if (file.exists(out_file)) {
      cat("  on disk:", basename(out_file), "\n")
      next
    }
    reps_b <- min(block_sz, n_sim - (b - 1L) * block_sz)

    set.seed(cl_spec$seed + b)
    res <- do.call(run_coverage_sim,
                   c(list(n_sim = reps_b, n_cores = n_cores), cl_spec$args))

    res$sweep <- cl_spec$sweep
    res$key   <- cl_spec$key
    res$block <- b
    for (nm in names(cl_spec$extra)) res[[nm]] <- cl_spec$extra[[nm]]

    # Write beside the target and rename: a task killed mid-write then leaves
    # no half-written file for the merge to read.
    tmp <- paste0(out_file, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, out_file)
    cat("  wrote:", basename(out_file), "\n")
  }
}

cat("task", array_id, "done\n")
