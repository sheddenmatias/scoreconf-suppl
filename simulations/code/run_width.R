# Run the width cells assigned to one SLURM array task; see run_cell.R for
# the layout conventions, which this follows: strided cells, blocks written
# as they finish, blocks on disk skipped, block seeds independent of the
# worker count.
#
# Environment:
#   SLURM_ARRAY_TASK_ID, SLURM_ARRAY_TASK_COUNT  set by sbatch
#   SLURM_CPUS_PER_TASK                          workers for the inner loop
#   RECONF_DIR                                   where the scripts are
#   RECONF_BLOCK                                 replications per block (250)
#   RECONF_WIDTHDIR                              block directory
#   RECONF_NSIM                                  replications per cell (1000)

root <- normalizePath(Sys.getenv("RECONF_DIR", unset = getwd()))
setwd(root)

locate <- function(name) {
  for (p in c(file.path(root, name),
              file.path(root, "R", name),
              file.path(root, "simulations", "code", name)))
    if (file.exists(p)) return(p)
  stop("cannot find ", name, " under ", root)
}
here <- function(...) locate(basename(file.path(...)))

source(locate("coverage_engine.R"))
source(locate("width_engine.R"))

env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(v) || v < 1L) default else v
}

array_id <- env_int("SLURM_ARRAY_TASK_ID", 1L)
n_arrays <- env_int("SLURM_ARRAY_TASK_COUNT", 1L)
n_cores  <- env_int("SLURM_CPUS_PER_TASK",
                    env_int("RECONF_CORES", max(1L, parallel::detectCores() - 1L)))
block_sz <- env_int("RECONF_BLOCK", 250L)
n_sim    <- env_int("RECONF_NSIM", 1000L)

width_dir <- Sys.getenv("RECONF_WIDTHDIR",
                        unset = file.path(root, "results", "width_blocks"))
dir.create(width_dir, recursive = TRUE, showWarnings = FALSE)

cells <- width_cells(n_sim)
pos   <- seq_along(cells)
mine  <- pos[(pos - array_id) %% n_arrays == 0]

cat(sprintf("task %d of %d: %d cells, %d cores, blocks of %d\n",
            array_id, n_arrays, length(mine), n_cores, block_sz))

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
  cell     <- cells[[ci]]
  n_blocks <- ceiling(cell$n_sim / block_sz)
  for (b in seq_len(n_blocks)) {
    out_file <- file.path(width_dir, sprintf("%s__b%02d.Rds", cell$key, b))
    if (file.exists(out_file)) {
      cat("  on disk:", basename(out_file), "\n")
      next
    }
    reps_b <- min(block_sz, cell$n_sim - (b - 1L) * block_sz)
    res <- run_width_block(cell, reps_b, cell$seed + b)
    if (is.null(res)) next
    res$block <- b
    tmp <- paste0(out_file, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, out_file)
    cat("  wrote:", basename(out_file), "\n")
  }
}

cat("task", array_id, "done\n")
