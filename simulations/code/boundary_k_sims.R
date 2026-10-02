library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
library(here)

# Simulation for Figure 3. The four statistics are computed by the
# self-contained engine in simulations/code/boundary_k_model.R.
#
# Environment (all optional; run_local.R sets them):
#   RECONF_NSIM    replications per setting (default 10000, as in the paper)
#   RECONF_CORES   parallel workers (default: all cores but one)
#   RECONF_RESDIR  where boundary_k_dat_both.Rds is written (default results/)
# The designs are always read from results/boundary_k_designs.Rds.

source(here("simulations", "code", "boundary_k_model.R"))

env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(v) || v < 1L) default else v
}

## Coverage simulation for one design
#
# Arguments:
#   n_sim  : number of replications
#   design : list with A (design matrix), u (stratum weights), R (target
#            information correlation), k
#   dfper  : degrees of freedom per unit weight; N_m = round(dfper * u_m)
#
# Returns a data frame with cols: k, n, stat, prop (coverage), se, reps
run_boundary_sim <- function(n_sim, design, dfper, eps = 0,
                             n_cores = max(1L, parallel::detectCores() - 1L)) {
  A <- design$A
  k <- design$k
  N <- round(dfper * design$u)

  # Self-check: the geometry actually simulated (from the integer N_m) must
  # match design$R, the geometry the cone prediction uses. Bounded weights keep
  # them equal; stop rather than report a discretization-distorted result.
  lik  <- bnd_mk_lik(A, N)
  keep <- 1:(k + 1)
  Ich  <- lik$info(c(1, rep(0, k), 1))
  Rchk <- cov2cor(Ich[keep, keep] - tcrossprod(Ich[keep, k + 2]) / Ich[k + 2, k + 2])
  if (max(abs(Rchk - design$R)) > 0.01)
    stop(sprintf("k=%d: realized geometry != design$R (max dev %.3f)",
                 k, max(abs(Rchk - design$R))))

  # Nuisances at eps, not on the boundary; the check above stays at zero
  # because design$R and the cone prediction are on-boundary quantities.
  psi_true <- c(1, rep(eps, k), 1)

  one_rep <- function(...) {
    SS <- bnd_sim(A, N, psi_true = psi_true)
    tryCatch(
      bnd_all_stats(A, N, SS, k),
      error = function(e) c(WLD = NA_real_, LRT = NA_real_,
                            SCR_C = NA_real_, SCR_U = NA_real_)
    )
  }

  cat(sprintf("[%s] Starting: k=%d n=%d  (%d reps, %d cores)\n",
              format(Sys.time(), "%H:%M:%S"), k, sum(N), n_sim, n_cores))
  t0 <- proc.time()["elapsed"]

  # Reproducible parallel RNG on a reused PSOCK cluster, exactly as in
  # coverage_engine.R (doRNG: one stream per iteration, backend independent;
  # PSOCK avoids fork, safe with Accelerate BLAS and portable to Windows).
  cl <- getOption("reconf.cluster")
  own_cluster <- is.null(cl)
  if (own_cluster) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    invisible(parallel::clusterEvalQ(cl, Sys.setenv(VECLIB_MAXIMUM_THREADS = "1",
                                                    OMP_NUM_THREADS = "1")))
    parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
  }
  results <- foreach::foreach(i = seq_len(n_sim)) %dorng% one_rep(i)
  if (own_cluster) parallel::stopCluster(cl)
  mat <- do.call(cbind, results)   # 4 x n_sim

  cat(sprintf("  Done. Total: %s\n\n", format(round(proc.time()["elapsed"] - t0))))

  q_crit <- qchisq(0.95, 1)
  prop   <- rowMeans(mat <= q_crit, na.rm = TRUE)
  reps   <- rowSums(!is.na(mat))
  se     <- sqrt(prop * (1 - prop) / reps)
  data.frame(k = design$k, n = sum(N), stat = names(prop),
             prop = prop, se = se, reps = reps, row.names = NULL)
}

###############
# SIMULATIONS #
###############
ks          <- 1:5
dfper       <- 20
n_sim_final <- env_int("RECONF_NSIM", 10000L)
n_design    <- 64      # random starts for the design search
eps_nuis    <- 0.001   # nuisances just off their boundary
kappa_fixed <- 0.99    # multiple correlation, held across k
n_target    <- 800     # sample size, held across k and across senses
senses      <- c(min = "min", max = "max")

.n_cores <- env_int("RECONF_CORES", max(1L, parallel::detectCores() - 1L))
.cl <- parallel::makeCluster(.n_cores)
doParallel::registerDoParallel(.cl)
invisible(parallel::clusterEvalQ(.cl, Sys.setenv(VECLIB_MAXIMUM_THREADS = "1",
                                                 OMP_NUM_THREADS = "1")))
parallel::clusterExport(.cl, ls(globalenv()), envir = globalenv())
options(reconf.cluster = .cl)

# Designs (cached). For each k and each sense, search loadings/weights to
# realize the extreme geometry with multiple correlation kappa_fixed, then
# rescale the weights to a common n. Seeded independently of later phases so a
# cached run and a regenerated run give the same coverage results.
design_file <- here("results", "boundary_k_designs.Rds")
if (file.exists(design_file)) {
  designs <- readRDS(design_file)
} else {
  designs <- lapply(senses, function(sn) {
    set.seed(11)
    Rtargets <- lapply(ks, function(k) bnd_extreme_R(k, kappa_fixed, sense = sn))
    names(Rtargets) <- as.character(ks)
    set.seed(7)
    out <- lapply(ks, function(k) {
      Rt   <- Rtargets[[as.character(k)]]
      cand <- foreach::foreach(s = seq_len(n_design)) %dorng% bnd_design_candidate(k, Rt)
      best <- cand[[which.min(sapply(cand, `[[`, "loss"))]]
      bnd_rescale_n(bnd_finalize_design(best$z, k, best$M), dfper, n_target)
    })
    names(out) <- as.character(ks)
    out
  })
  saveRDS(designs, design_file)
}

run_sense <- function(sn) {
  ds <- designs[[sn]]
  set.seed(3)
  cone <- vapply(ks, function(k) 1 - bnd_cone_noncov(ds[[as.character(k)]]$R),
                 numeric(1))
  cat(sprintf("=== boundary-nuisance coverage, sense = %s, r = 3..7 ===\n", sn))
  set.seed(101)
  stats <- do.call(rbind, lapply(ks, function(k)
    run_boundary_sim(n_sim = n_sim_final, design = ds[[as.character(k)]],
                     dfper = dfper, eps = eps_nuis, n_cores = .n_cores)))
  # Carry the cone prediction as its own series (se = 0) for the figure.
  cone_df <- data.frame(k = ks, n = NA_integer_, stat = "CONE", prop = cone,
                        se = 0, reps = NA_integer_, row.names = NULL)
  cbind(rbind(stats, cone_df), sense = sn)
}

boundary_dat <- do.call(rbind, lapply(senses, run_sense))

res_dir <- Sys.getenv("RECONF_RESDIR", unset = here("results"))
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(boundary_dat, file = file.path(res_dir, "boundary_k_dat_both.Rds"))

parallel::stopCluster(.cl)
options(reconf.cluster = NULL)
