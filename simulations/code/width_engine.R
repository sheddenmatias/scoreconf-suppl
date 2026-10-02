# Interval endpoints for the width comparison. Requires coverage_engine.R to
# be sourced first, for the simulators, the constrained optimizers and the
# index mapping.
#
# Five constructions of a confidence interval for the random-effect
# covariance: Wald from the efficient observed information (absolute value),
# likelihood ratio and score on the constrained set (both from one outward
# search over constrained null fits), and the score on the extended set with
# full and one-step nuisance estimation via reconf::ci_lmer.

WIDTH_Q_CRIT <- qchisq(0.95, 1)

# The cells used in the paper, with the keys and seeds of the original run.
# The covariance psi2 is the parameter of interest throughout, so every bound is sought over the whole real line.
#
#   snr_20 ... snr_160  Figure 2: psi = (2, 0.4, 2, 1) with 20, 40, 80
#                       or 160 groups of 10 observations.
#   rho_0.999           the Section 4 statement that the constrained score
#                       interval is empty in 54% of replications at
#                       correlation 0.999: psi = (2.5, 2.5 * 0.999, 2.5, 1),
#                       20 groups of 10.
width_cells <- function(n_sim = 1000L) {
  cells <- list(
    list(key = "rho_0.999", psi = c(2.5, 2.5 * 0.999, 2.5, 1),
         n1 = 20L, n2 = 10L, idx = 2L, seed = 1320007L, n_sim = n_sim))
  snr_n1 <- c(20L, 40L, 80L, 160L)
  for (i in seq_along(snr_n1)) {
    cells[[length(cells) + 1L]] <-
      list(key = sprintf("snr_%d", snr_n1[i]), psi = c(2, 0.4, 2, 1),
           n1 = snr_n1[i], n2 = 10L, idx = 2L, seed = 1340000L + i,
           n_sim = n_sim)
  }
  names(cells) <- vapply(cells, `[[`, character(1), "key")
  cells
}

# Endpoints of the likelihood-ratio and constrained-score intervals: from the
# constrained full maximum outward in steps of hw/4, one constrained null fit
# per step warm started at the previous solution, each statistic's bound by
# linear interpolation at its first crossing of WIDTH_Q_CRIT.
constrained_ci <- function(fit, fit_idx, hw) {
  Y     <- getME(fit, "y")
  X     <- getME(fit, "X")
  Z     <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)

  psi_hat     <- reconf:::get_psi_hat_lmer(fit)
  psd_block   <- corr_block_idx(fit)
  if (is.null(psd_block)) stop("the fit has no correlated random-effect block")
  free_idx    <- seq_along(psi_hat)[-fit_idx]

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

  # Constrained full maximum, as in all_stats: lme4's estimate plus jittered
  # restarts, so the likelihood ratio is non-negative by construction.
  full_starts <- c(list(psi_hat),
                   lapply(1:3, function(j)
                     psi_hat * runif(length(psi_hat), 0.5, 2)))
  fit_full <- best_of(full_starts, seq_along(psi_hat))
  if (is.null(fit_full)) return(NULL)
  ll_full <- fit_full$value

  stat_at <- function(v, warm) {
    s0 <- warm
    s0[fit_idx] <- v
    fc <- cfit(s0, free_idx)
    if (is.null(fc)) return(NULL)
    list(lrt = max(0, 2 * (ll_full - fc$value)),
         scr = tryCatch(
           reconf:::score_stat(theta = fc$arg, test_idx = fit_idx,
                               Y = Y, X = X, Z = Z, Hlist = Hlist,
                               expected = TRUE, REML = TRUE,
                               signed = FALSE, efficient = TRUE),
           error = function(e) NA_real_),
         arg = fc$arg)
  }

  center <- fit_full$arg[fit_idx]
  st0    <- stat_at(center, fit_full$arg)
  if (is.null(st0)) return(NULL)

  search <- function(dir) {
    bounds <- c(lrt = NA_real_, scr = NA_real_)
    # A statistic already past the critical value at the center has no
    # crossing to find.
    for (s in names(bounds))
      if (!is.na(st0[[s]]) && st0[[s]] >= WIDTH_Q_CRIT) bounds[[s]] <- center
    step   <- hw / 4
    v      <- center
    warm   <- fit_full$arg
    prev   <- st0
    prev_v <- v
    for (k in 1:400) {
      if (!anyNA(bounds)) break
      v  <- v + dir * step
      st <- stat_at(v, warm)
      if (is.null(st)) break
      warm <- st$arg
      for (s in names(bounds)) {
        if (is.na(bounds[[s]]) && !is.na(st[[s]]) && !is.na(prev[[s]]) &&
            st[[s]] >= WIDTH_Q_CRIT)
          bounds[[s]] <- prev_v + (v - prev_v) *
            (WIDTH_Q_CRIT - prev[[s]]) / (st[[s]] - prev[[s]])
      }
      prev   <- st
      prev_v <- v
      if (k %% 50 == 0) step <- 2 * step
    }
    bounds
  }

  up <- search(+1)
  lo <- search(-1)
  list(lrt = c(lo[["lrt"]], up[["lrt"]]),
       scr = c(lo[["scr"]], up[["scr"]]))
}

# One replication: simulate, compute the five intervals, return one row per
# statistic with endpoints.
width_one_rep <- function(psi, n1, n2, idx) {
  fit     <- suppressWarnings(sim_corr(psi, n1, n2, p = 1))
  fit_idx <- grp_psi_idx(fit, idx, "corr")

  Y     <- getME(fit, "y")
  X     <- getME(fit, "X")
  Z     <- getME(fit, "Z")
  Hlist <- reconf:::get_Hlist_lmer(fit)
  psi_hat <- reconf:::get_psi_hat_lmer(fit)

  # Wald from the efficient observed information, absolute value as in
  # all_stats.
  ll_hat <- tryCatch(
    reconf:::loglikelihood(psi = psi_hat, Y = Y, X = X, Z = Z, Hlist = Hlist,
                           REML = TRUE, get_val = FALSE, get_score = FALSE,
                           get_inf = TRUE, expected = FALSE),
    error = function(e) NULL)
  wald <- c(NA_real_, NA_real_)
  hw   <- NA_real_
  if (!is.null(ll_hat)) {
    I    <- ll_hat$inf_mat
    ei   <- drop(I[fit_idx, fit_idx] -
                   I[fit_idx, -fit_idx, drop = FALSE] %*%
                   solve(I[-fit_idx, -fit_idx, drop = FALSE],
                         I[-fit_idx, fit_idx, drop = FALSE]))
    hw   <- sqrt(WIDTH_Q_CRIT / abs(ei))
    wald <- psi_hat[fit_idx] + c(-1, 1) * hw
  }

  # Score on the extended set, full and one-step nuisance estimation. The
  # search radius is set wide enough to bracket the bounds in every cell.
  ci_row <- function(...) tryCatch({
    ci <- reconf::ci_lmer(fit, test_idx = fit_idx, level = 0.95,
                          expected = TRUE, num_points = 20000L, ...)
    c(ci[1, "lower"], ci[1, "upper"])
  }, error = function(e) c(NA_real_, NA_real_))
  scr_u <- ci_row(onestep = FALSE)
  scr_1 <- ci_row(onestep = TRUE)

  # Likelihood ratio and score on the constrained set.
  cci <- if (is.na(hw)) NULL else constrained_ci(fit, fit_idx, hw)
  lrt   <- if (is.null(cci)) c(NA_real_, NA_real_) else cci$lrt
  scr_c <- if (is.null(cci)) c(NA_real_, NA_real_) else cci$scr

  out <- rbind(WLD = wald, LRT = lrt, SCR_C = scr_c, SCR_U = scr_u,
               SCR_1 = scr_1)
  data.frame(stat = rownames(out), lower = out[, 1], upper = out[, 2],
             row.names = NULL)
}

# One block of replications for one cell, on the registered parallel backend.
run_width_block <- function(cell, n_sim, seed) {
  cat(sprintf("[%s] Starting: %s psi=(%s)  (%d reps)\n",
              format(Sys.time(), "%H:%M:%S"), cell$key,
              paste(round(cell$psi, 4), collapse = ","), n_sim))
  t0 <- proc.time()["elapsed"]
  set.seed(seed)
  res <- foreach::foreach(
    i = seq_len(n_sim),
    .packages = c("Matrix", "lme4", "nloptr", "alabama", "reconf")
  ) %dorng% {
    out <- tryCatch(width_one_rep(cell$psi, cell$n1, cell$n2, cell$idx),
                    error = function(e) NULL)
    if (!is.null(out)) out$rep <- i
    out
  }
  cat(sprintf("  Done. Total: %s\n", format(round(proc.time()["elapsed"] - t0))))
  out <- do.call(rbind, res)
  if (is.null(out)) return(NULL)
  out$key  <- cell$key
  out$true <- cell$psi[cell$idx]
  out
}
