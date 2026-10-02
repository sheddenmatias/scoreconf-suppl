# Shared engine for the lme4-based simulations (coverage for Figure 1 and the width
# comparison for Figure 2): the constrained fits, the four statistics, and one replication.
# Defines: maximize_loglik_3(), all_stats(), run_coverage_sim()

library(Matrix)
library(lme4)
library(alabama)
library(parallel)
library(nloptr)
library(here)
library(foreach)
library(doParallel)
library(doRNG)

# Simulations

source(here("simulations", "code", "simulators.R"))

# Constrained optimization for the CORR model 
maximize_loglik_3 <- function(start_val, opt_idx, Y, X, Z, Hlist, expected = TRUE,
                              REML = TRUE, precomp = NULL, blk = 1:3) {
  n <- length(Y)
  r <- length(Hlist) + 1
  p <- ncol(X)
  expected_length <- if(REML) r else p + r
  
  # reconf >= 0.1.1: get_precomp also caches H and the factor used by the
  # feasibility check, so repeated likelihood evaluations stay fast
  precomp <- reconf:::get_precomp(Y = Y, X = X, Z = Z, REML = TRUE,
                                  Hlist = Hlist)
  
  # objective function
  obj_fun <- function(x) {
    psi <- start_val
    psi[opt_idx] <- x
    # reconf >= 0.1.1: loglikelihood() builds Psi, checks feasibility (value
    # -Inf when Sigma is not positive definite), and factorizes internally
    ll_things <- reconf:::loglikelihood(psi = psi, Y = Y, X = X, Z = Z,
                                        Hlist = Hlist, REML = TRUE,
                                        precomp = precomp, check = FALSE)
    list("value" = -ll_things$value, "gradient" = -ll_things$score[opt_idx],
         "hessian" = as.matrix(ll_things$inf_mat[opt_idx, opt_idx]))
  }
  
  # separate function and its gradient
  f_val <- function(x) {
    val <- obj_fun(x)$value
    if (!is.finite(val)) return(1e10)
    return(val)
  }
  f_grad <- function(x) {
    g <- obj_fun(x)$gradient
    if (any(!is.finite(g))) rep(0, length(x)) else g
  }
  
  # param of interest
  idx <- seq_along(start_val)[-opt_idx]
  covar <- start_val[idx]
  
  # Constraints for the correlated random-effects model, written in terms of the
  # FULL parameter vector so that any subset may be held fixed. With
  # psi = (var intercept, covariance, var slope, error variance) the parameter
  # set is
  #     psi1 >= 0,  psi3 >= 0,  psi1 psi3 - psi2^2 >= 0,  psi4 >= 0,
  # and psi2 is unconstrained in sign, a covariance having no boundary of its
  # own.
  psi_full <- function(x) { psi <- start_val; psi[opt_idx] <- x; psi }

  # psi[1:3] are the correlated block (var, covariance, var); psi[4:r] are
  # variances of INDEPENDENT random effects and the error variance, each of
  # which needs only non-negativity.
  # blk gives the psi positions of the correlated block (variance, covariance,
  # variance); every remaining parameter is a variance of an independent random
  # effect, or the error variance, and needs only non-negativity. lme4 orders
  # grouping factors by number of levels, so blk is not always 1:3.
  oth <- setdiff(seq_len(r), blk)

  hin <- function(x) {
    psi <- psi_full(x)
    c(psi[blk[1]] * psi[blk[3]] - psi[blk[2]]^2, psi[blk[1]], psi[blk[3]],
      psi[oth])
  }
  hin.jac <- function(x) {
    psi <- psi_full(x)
    J <- matrix(0, 3 + length(oth), r)
    J[1, blk] <- c(psi[blk[3]], -2 * psi[blk[2]], psi[blk[1]])
    J[2, blk[1]] <- 1
    J[3, blk[3]] <- 1
    for (j in seq_along(oth)) J[3 + j, oth[j]] <- 1
    J[, opt_idx, drop = FALSE]
  }
  
  hin2 <- function(x) -hin(x)
  hin.jac2 <- function(x) -hin.jac(x)
  
  # optional
  # helpful if information matrix is nearly computationally singular
  find_feasible_start <- function(x_init, covar, min_var = 1e-8) {
    
    if (all(hin(x_init) >= 0)) return(x_init)
    
    message("Initial value infeasible, finding feasible start...")
    
    x_safe <- pmax(x_init, min_var)
    
    proj <- alabama::auglag(
      par           = x_safe,
      fn            = function(x) sum((x - x_safe)^2),  # minimize distance from x_init
      gr            = function(x) 2*(x - x_safe),
      hin           = hin,
      hin.jac       = hin.jac,
      control.outer = list(trace = FALSE, eps = 1e-4),
      control.optim = list(trace = 0)
    )
    
    if (!all(hin(proj$par) >= 0)) {
      warning("Could not find feasible start, using fallback")
      return(c(2*abs(covar), 2*abs(covar), x_init[3]))  # safe default
    }
    
    proj$par
  }
  # can use feasible_start in what follows
  # feasible_start <- find_feasible_start(start_val[opt_idx], covar)
  
  fit <- slsqp(x0 = start_val[opt_idx],
               fn = f_val,
               hin = hin2,
               hinjac = hin.jac2,
               deprecatedBehavior = FALSE)
  fit_par <- fit$par
  
  if (fit$convergence < 0) {
    warning("Optimization did not converge. Results may be unreliable. ",
            "Iterations: ", fit$counts)
    warning(fit$message)
  }
  # Return results
  start_val[opt_idx] <- fit_par
  names(start_val) <- if(REML) paste0("psi", 1:r) else c(paste0("b", 1:p), paste0("psi", 1:r))
  list("arg" = start_val, "value" = -fit$value, "conv" = fit$convergence,
       "iter" = fit$counts)
}

## All test statistics in one pass
#
#   - Y, X, Z, Hlist extracted once
#   - psi_hat and ll_hat (loglikelihood at MLE, value + info) computed once
#     - used by Wald (info matrix) and LR (value)
#   - constrained null fit (nuisance optimized, idx fixed at null value) computed once
#     - used by LR (value) and SCR_C (arg)
#   - unconstrained null fit (reconf:::maximize_loglik, trust region) computed once
#     - used by SCR_U
#
# The constrained fits use maximize_loglik_3, which keeps the random-effect
# covariance matrix positive semidefinite.

# Feasible starting point: projects a parameter vector into the proper parameter
# set, moving only the coordinates that are free to move. Near the boundary the 
# natural starting value for the null fit, the fitted estimate with the interest
# parameter replaced by its null value, is frequently infeasible.
make_feasible <- function(psi, opt_idx, blk, eps = 1e-6) {
  # variances must be non-negative; a covariance has no such constraint
  vars <- setdiff(opt_idx, if (is.null(blk)) integer(0) else blk[2])
  psi[vars] <- pmax(psi[vars], eps)
  if (!is.null(blk)) {
    b1 <- blk[1]; b2 <- blk[2]; b3 <- blk[3]
    need <- psi[b2]^2
    if (psi[b1] * psi[b3] < need) {
      if (b1 %in% opt_idx && b3 %in% opt_idx) {
        psi[b1] <- max(psi[b1], sqrt(need) + eps)
        psi[b3] <- max(psi[b3], need / psi[b1] + eps)
      } else if (b1 %in% opt_idx) {
        psi[b1] <- need / max(psi[b3], eps) + eps
      } else if (b3 %in% opt_idx) {
        psi[b3] <- need / max(psi[b1], eps) + eps
      } else if (b2 %in% opt_idx) {
        psi[b2] <- sign(psi[b2]) * max(sqrt(max(psi[b1] * psi[b3], 0)) - eps, 0)
      }
    }
  }
  psi
}

all_stats <- function(fit_lme4, psi_start, idx, type) {
  Y     <- getME(fit_lme4, "y")
  X     <- getME(fit_lme4, "X")
  Z     <- getME(fit_lme4, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit_lme4)

  # Usual MLE (constrained), its restricted log-likelihood for the LR, and the
  # OBSERVED information for Wald -- what the common Wald statistic uses.
  psi_hat <- reconf:::get_psi_hat_lmer(fit_lme4)
  ll_hat  <- reconf:::loglikelihood(psi      = psi_hat,
                                    Y        = Y, X = X, Z = Z, Hlist = Hlist,
                                    REML     = TRUE,
                                    get_val  = TRUE, get_score = FALSE,
                                    get_inf  = TRUE, get_beta  = FALSE,
                                    expected = FALSE)

  # constrained null fit: nuisance params optimized, idx held at null (shared by LR and SCR_C)
  psd_block   <- corr_block_idx(fit_lme4)
  if (is.null(psd_block)) stop("the fit has no correlated random-effect block")
  free_idx    <- seq_along(psi_start)[-idx]

  # One constrained maximization from a given start, over a given free set.
  cfit <- function(s0, oidx) {
    s0 <- make_feasible(s0, oidx, psd_block)
    tryCatch(
      maximize_loglik_3(start_val = s0, opt_idx = oidx,
                        Y = Y, X = X, Z = Z, Hlist = Hlist, REML = TRUE,
                        blk = psd_block),
      error = function(e) NULL)
  }
  best_of <- function(starts, oidx) {
    best <- NULL
    for (s0 in starts) {
      o <- cfit(s0, oidx)
      if (!is.null(o) && (is.null(best) || o$value > best$value)) best <- o
    }
    best
  }

  # Null fit: several feasible starts, keep the best.
  null_starts <- c(list(psi_start),
                   lapply(1:3, function(j) {
                     s <- psi_start
                     s[free_idx] <- s[free_idx] * runif(length(free_idx), 0.5, 2)
                     s
                   }))
  fit_null_c <- best_of(null_starts, free_idx)

  # Full fit: lme4's estimate need not maximize over the proper parameter set
  # near the boundary. Re-maximize from the null fit over the SAME set and keep
  # the better of the two. Because the null is a restriction of the full model,
  # this makes the likelihood ratio statistic non-negative by construction, as
  # in boundary_k_model.R.
  fit_full_c <- best_of(list(fit_null_c$arg, psi_hat), seq_along(psi_start))
  ll_full    <- max(ll_hat$value,
                    if (is.null(fit_full_c)) -Inf else fit_full_c$value)

  # unconstrained null fit (only for SCR_U). psi_start holds the tested
  # component at its null value alongside fitted nuisances, a combination that
  # can lie outside the extended set; trust() then rejects the start, so it is
  # made feasible first. A failure here is confined to SCR_U.
  ufit <- function(s0) tryCatch(
    reconf:::maximize_loglik(start_val = s0, opt_idx = free_idx,
                             Y = Y, X = X, Z = Z, Hlist = Hlist,
                             expected  = TRUE, REML = TRUE),
    error = function(e) NULL)
  # Same multi-start as the constrained fit.
  u_starts <- c(list(make_feasible(psi_start, free_idx, psd_block)),
                lapply(1:3, function(j) {
                  s <- psi_start
                  s[free_idx] <- s[free_idx] * runif(length(free_idx), 0.5, 2)
                  make_feasible(s, free_idx, psd_block)
                }))
  fit_null_u <- NULL
  for (s0 in u_starts) {
    o <- ufit(s0)
    if (!is.null(o) && (is.null(fit_null_u) || o$value > fit_null_u$value))
      fit_null_u <- o
  }

  # --- Wald ---
  inf_mat <- ll_hat$inf_mat[idx, idx, drop = FALSE]
  A_nt    <- ll_hat$inf_mat[-idx,  idx, drop = FALSE]
  I_nn    <- ll_hat$inf_mat[-idx, -idx, drop = FALSE]
  ei      <- inf_mat - crossprod(A_nt, solve(I_nn, A_nt))
  # The interest parameter is scalar, so the efficient information is a number.
  # Its absolute value is used, since the observed information can be negative
  # near the boundary.
  wld     <- (psi_hat[idx] - psi_start[idx])^2 * abs(ei)

  # --- LR ---
  lrt <- max(0, 2 * (ll_full - fit_null_c$value))

  # --- Score (constrained null) ---
  scr_c <- reconf:::score_stat(theta     = fit_null_c$arg,
                                test_idx  = idx,
                                Y = Y, X = X, Z = Z, Hlist = Hlist,
                                expected  = TRUE, REML = TRUE,
                                signed    = FALSE, efficient = TRUE)

  # --- Score (unconstrained null) ---
  scr_u <- if (is.null(fit_null_u)) NA_real_ else
    reconf:::score_stat(theta     = fit_null_u$arg,
                        test_idx  = idx,
                        Y = Y, X = X, Z = Z, Hlist = Hlist,
                        expected  = TRUE, REML = TRUE,
                        signed    = FALSE, efficient = TRUE)

  # NEGINF is not a statistic; it records that the observed information was
  # negative, so that run_coverage_sim can report the rate. Dropped there.
  c(WLD = wld, LRT = lrt, SCR_C = scr_c, SCR_U = scr_u,
    NEGINF = as.numeric(ei < 0))
}


## Coverage Probability Simulation
#
# Arguments:
#   n_sim   : number of replications
#   psi     : true parameter vector
#   n1, n2  : group sizes
#   n_per   : unused (kept for the argument lists in cells.R)
#   type    : "corr", the correlated random intercept and slope model
#   idx     : index of the parameter under test (default 1)
#   p       : number of non-intercept covariates
#   n_cores : number of cores for parallelization (default: detectCores() - 1)
#
# Returns data frame with cols:
#   type, n1, n2, p, psi1, ..., stat, prop, se, reps

run_coverage_sim <- function(n_sim, psi, n1, n2, n_per = NA, type,
                             idx = 1, p = 1,
                             n_cores = max(1L, parallel::detectCores() - 1L)) {
  r  <- length(psi)
  df <- 1

  one_rep <- function(...) {
    fit <- suppressWarnings(switch(type,
      corr = sim_corr(psi, n1, n2, p = p),
      stop("unknown simulation type: ", type)
    ))
    psi_start           <- reconf:::get_psi_hat_lmer(fit)
    fit_idx             <- grp_psi_idx(fit, idx, type)
    psi_start[fit_idx]  <- psi[idx]

    tryCatch(
      all_stats(fit, psi_start, fit_idx, type),
      error = function(e) c(WLD = NA_real_, LRT = NA_real_,
                             SCR_C = NA_real_, SCR_U = NA_real_,
                             NEGINF = NA_real_)
    )
  }

  psi_label <- paste0("psi=(", paste(round(psi, 3), collapse = ","), ")")
  cat(sprintf("[%s] Starting: type=%s n1=%d n2=%d %s  (%d reps, %d cores)\n",
              format(Sys.time(), "%H:%M:%S"), type, n1, n2, psi_label, n_sim, n_cores))
  t0 <- proc.time()["elapsed"]

  # Reproducible parallel RNG (doRNG: one stream per iteration, independent
  # of worker count and backend) on a PSOCK cluster (no fork: safe with
  # Accelerate BLAS on macOS, portable to Windows). Reuses the cluster that
  # run_cell.R sets up; creates one per call otherwise. NOTE: assumes the script is sourced into the global
  # environment (as Rscript does) so helpers reach the workers.
  cl <- getOption("reconf.cluster")
  own_cluster <- is.null(cl)
  if (own_cluster) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    invisible(parallel::clusterEvalQ(cl, Sys.setenv(VECLIB_MAXIMUM_THREADS = "1",
                                                    OMP_NUM_THREADS = "1")))
    parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
  }
  results <- foreach::foreach(
    i = seq_len(n_sim),
    .packages = c("Matrix", "lme4", "nloptr", "alabama", "reconf")
  ) %dorng% one_rep(i)
  if (own_cluster) parallel::stopCluster(cl)
  mat     <- do.call(cbind, results)
  n_neg_inf <- sum(mat["NEGINF", ] > 0, na.rm = TRUE)
  mat     <- mat[setdiff(rownames(mat), "NEGINF"), , drop = FALSE]

  cat(sprintf("  Done. Total: %s\n\n",
              format(round(proc.time()["elapsed"] - t0))))

  q_crit <- qchisq(0.95, df)
  prop   <- rowMeans(mat <= q_crit, na.rm = TRUE)
  reps   <- rowSums(!is.na(mat))
  se     <- sqrt(prop * (1 - prop) / reps)

  if (n_neg_inf > 0)
    cat(sprintf("  Wald: observed information negative in %d of %d (%.1f%%)\n",
                n_neg_inf, reps[["WLD"]], 100 * n_neg_inf / reps[["WLD"]]))

  psi_cols <- setNames(as.list(psi), paste0("psi", seq_len(r)))
  cbind(
    data.frame(type = type, n1 = n1, n2 = n2, p = p),
    as.data.frame(psi_cols),
    data.frame(stat = names(prop), prop = prop, se = se, reps = reps,
               row.names = NULL)
  )
}

