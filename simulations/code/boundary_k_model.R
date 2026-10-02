# Model, design search, and statistics for the coverage vs. 
# number-of-nuisance params simulation (Figure 3)
#
# Sourced by (simulations/code/boundary_k_sims.R) which runs the simulation.
# The model is a variance components model: M independent strata, stratum m has 
# N_m iid N(0, sigma_m) observations
# with sigma_m = sum_j psi_j h_jm + psi_res and known nonnegative loadings h_jm.
# Reducing each stratum to its sum of squares SS_m ~ sigma_m * chisq(N_m) makes
# large samples cheap. The loadings and stratum weights are searched to realize
# a target information geometry supplied by the user. This is NOT an
# lme4 random-effects design, so the
# statistics are computed here rather than through reconf; the definitions of
# WLD, LRT, SCR_C, SCR_U match all_stats() in simulations/code/coverage_engine.R.

## ---------------------------------------------------------------------------
## Information geometry and design search
## ---------------------------------------------------------------------------

# Partialled information correlation among (interest, nuisances) at the truth
# (nuisances = 0), after profiling out the interior residual variance.
bnd_Rstar <- function(hint, u, Hnu) {
  sig <- hint + 1
  H   <- cbind(hint, Hnu, 1)
  d   <- ncol(Hnu) + 2
  I   <- matrix(0, d, d)
  for (i in 1:d) for (j in i:d)
    I[i, j] <- I[j, i] <- 0.5 * sum(u * H[, i] * H[, j] / sig^2)
  kk    <- ncol(Hnu) + 1
  Istar <- I[1:kk, 1:kk] - tcrossprod(I[1:kk, d]) / I[d, d]
  list(R      = cov2cor(Istar),
       mineig = min(eigen(cov2cor(I), symmetric = TRUE, only.values = TRUE)$values))
}

# Exchangeable target: corr(interest, nuisance_j) = a, and tau the correlation
# between nuisance residuals after partialling out the interest score.
bnd_target_R <- function(k, a, tau) {
  d   <- k + 1
  R   <- diag(d); R[1, -1] <- a; R[-1, 1] <- a
  for (i in 1:k) for (j in 1:k) if (i != j) R[1 + i, 1 + j] <- a^2 + (1 - a^2) * tau
  R
}

# Multiple correlation between the interest score and the nuisance scores.
bnd_kappa <- function(R) {
  r <- R[1, -1, drop = FALSE]
  sqrt(drop(r %*% solve(R[-1, -1, drop = FALSE]) %*% t(r)))
}

# The tau that puts bnd_target_R(k, a, tau) at multiple correlation kappa.
# From kappa^2 = a^2 k / (1 + (k-1) rho) with rho = a^2 + (1 - a^2) tau.
bnd_tau_for <- function(a, k, kappa) {
  if (k == 1) return(0)
  (a^2 * k / kappa^2 - 1 - (k - 1) * a^2) / ((k - 1) * (1 - a^2))
}

# Over the pairwise correlation a, with tau set to hold the multiple
# correlation at kappa, the geometry minimizing ("min") or maximizing ("max")
# the cone model's predicted LRT coverage. At k = 1 the family has one member,
# since kappa = a there.
bnd_extreme_R <- function(k, kappa, sense = c("min", "max"),
                          a_grid = seq(0.50, 0.995, by = 0.005),
                          N = 2e5, min_eig = 1e-4) {
  sense <- match.arg(sense)
  if (k == 1) return(bnd_target_R(1, kappa, 0))
  best <- NULL
  for (a in a_grid) {
    if (a >= kappa) next
    tau <- bnd_tau_for(a, k, kappa)
    if (!is.finite(tau) || tau <= -1 / (k - 1) || tau >= 1) next
    R <- bnd_target_R(k, a, tau)
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < min_eig) next
    if (abs(bnd_kappa(R) - kappa) > 1e-6) next
    cov_k <- 1 - bnd_cone_noncov(R, N = N)
    better <- is.null(best) ||
      (sense == "min" && cov_k < best$cov) || (sense == "max" && cov_k > best$cov)
    if (better) best <- list(R = R, cov = cov_k)
  }
  if (is.null(best)) stop("no feasible geometry at k = ", k, ", kappa = ", kappa)
  best$R
}

# Scale the stratum weights so that n = dfper * sum(u) hits n_target. The
# target correlation matrix is invariant to a common scaling of u, so this
# changes the sample size and nothing else.
bnd_rescale_n <- function(design, dfper, n_target) {
  design$u <- design$u * n_target / (dfper * sum(design$u))
  design
}

# Bounded reparametrizations. Weights u stay in [0.25, 4] so that discretizing
# the degrees of freedom to integers (round(dfper * u)) does not distort the
# realized geometry.
bnd_ub <- function(x) 0.25 + 3.75 * plogis(x)
bnd_hb <- function(x) 0.02 * exp(log(50 / 0.02) * plogis(x))
bnd_nb <- function(x) 1e-3 * exp(log(5e4) * plogis(x))
bnd_decode <- function(z, M, k)
  list(hint = bnd_hb(z[1:M]), u = bnd_ub(z[M + 1:M]),
       Hnu = matrix(bnd_nb(z[2 * M + 1:(M * k)]), M, k))

# Distance from the achieved to the target geometry, penalizing near-singular
# information (keeps the design estimable).
bnd_design_loss <- function(z, k, M, Rtarget) {
  p <- bnd_decode(z, M, k)
  o <- bnd_Rstar(p$hint, p$u, p$Hnu)
  pen <- if (o$mineig < 1e-3) 30 * (1e-3 - o$mineig) / 1e-3 else 0
  sum((o$R - Rtarget)^2) + pen
}

# One random-start fit of the design to its target (RNG supplied by the caller).
bnd_design_candidate <- function(k, Rtarget, M = 3 * k + 5) {
  z0 <- rnorm(M * (k + 2), 0, 1.6)
  o  <- optim(z0, function(z) bnd_design_loss(z, k, M, Rtarget), method = "BFGS",
              control = list(maxit = 2500, reltol = 1e-12))
  list(loss = o$value, z = o$par, M = M)
}

# Assemble a design object (design matrix A and target correlation R) from an
# optimized parameter vector.
bnd_finalize_design <- function(z, k, M = 3 * k + 5) {
  p <- bnd_decode(z, M, k)
  list(k = k, M = M, u = p$u, A = cbind(p$hint, p$Hnu, 1),
       R = bnd_Rstar(p$hint, p$u, p$Hnu)$R)
}

## ---------------------------------------------------------------------------
## Cone-model prediction: the local-asymptotic law of the profile LRT when the
## nuisances sit on their boundary. Both fits maximize a quadratic over the
## nuisance orthant; the value is found by enumerating faces (the QP optimum is
## the largest primal-feasible stationary point since R is positive definite).
## ---------------------------------------------------------------------------

bnd_cone_value <- function(S, R, gamma_free) {
  N <- nrow(S); d <- ncol(R); lam <- 2:d; k <- length(lam); v <- numeric(N)
  for (m in 0:(2^k - 1)) {
    freeL <- lam[bitwAnd(m, 2^(0:(k - 1))) > 0]
    Fset  <- c(if (gamma_free) 1, freeL)
    if (!length(Fset)) next
    Tmat <- S[, Fset, drop = FALSE] %*% solve(R[Fset, Fset, drop = FALSE])
    li   <- which(Fset %in% lam)
    ok   <- if (length(li)) rowSums(Tmat[, li, drop = FALSE] < 0) == 0 else rep(TRUE, N)
    if (any(ok))
      v[ok] <- pmax(v[ok], 0.5 * rowSums(Tmat[ok, , drop = FALSE] * S[ok, Fset, drop = FALSE]))
  }
  v
}

# Predicted non-coverage of the nominal 95% LRT interval for a given geometry.
bnd_cone_noncov <- function(R, N = 4e5) {
  S <- matrix(rnorm(N * ncol(R)), N) %*% chol(R)
  mean(2 * (bnd_cone_value(S, R, TRUE) - bnd_cone_value(S, R, FALSE)) > qchisq(0.95, 1))
}

## ---------------------------------------------------------------------------
## Likelihood, score, expected information (sigma = A %*% psi)
## ---------------------------------------------------------------------------

bnd_mk_lik <- function(A, N) list(
  nll = function(psi, SS) { s <- as.vector(A %*% psi)
    if (any(s <= 0)) 1e12 else 0.5 * sum(N * log(s) + SS / s) },
  gr  = function(psi, SS) { s <- as.vector(A %*% psi)
    if (any(s <= 0)) rep(0, ncol(A)) else as.vector(crossprod(A, 0.5 * (N / s - SS / s^2))) },
  info = function(psi) { s <- as.vector(A %*% psi); 0.5 * t(A * (N / s^2)) %*% A }
)

# Robust box-constrained minimizer: many starts, KKT-checked, refit on failure.
# A lazy full-model fit would shrink the LRT and hide the effect, so convergence
# is verified rather than assumed.
bnd_fit_min <- function(nll, gr, lb, starts, ktol = 1e-4) {
  run <- function(st) tryCatch(
    optim(st, nll, gr, method = "L-BFGS-B", lower = lb,
          control = list(factr = 1e1, pgtol = 1e-10, maxit = 2000)),
    error = function(e) NULL)
  kkt_ok <- function(o) { g <- gr(o$par); at <- o$par <= lb + 1e-8
    all(abs(g[!at]) < ktol * (1 + abs(o$value))) &&
      all(g[at] > -ktol * (1 + abs(o$value))) }
  best <- NULL
  for (st in starts) { o <- run(st)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  tries <- 0
  while ((is.null(best) || !kkt_ok(best)) && tries < 3) {
    tries <- tries + 1
    for (j in 1:8) { st <- pmax(lb, abs(rnorm(length(lb), 1, 1)))
      o <- run(st); if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  }
  list(par = best$par, value = best$value, kkt = kkt_ok(best))
}

# Unrestricted (extended-set) nuisance maximizer for SCR_U: interest fixed, the
# k boundary components and residual free; self-barriered since nll -> Inf as
# any sigma_m -> 0.
bnd_fit_unrestricted <- function(A, N, SS, k, p1 = 1) {
  Ar <- A[, -1, drop = FALSE]; off <- p1 * A[, 1]
  nll <- function(e) { s <- off + as.vector(Ar %*% e)
    if (any(s <= 0)) 1e12 else 0.5 * sum(N * log(s) + SS / s) }
  gr  <- function(e) { s <- off + as.vector(Ar %*% e)
    if (any(s <= 0)) rep(0, k + 1) else as.vector(crossprod(Ar, 0.5 * (N / s - SS / s^2))) }
  starts <- c(list(c(rep(0, k), 1), rep(1, k + 1)),
              lapply(1:6, function(i) pmax(1e-3, rnorm(k + 1, 0.5, 0.5))))
  best <- NULL
  for (st in starts) {
    o1 <- tryCatch(optim(st, nll, method = "Nelder-Mead",
                         control = list(maxit = 3000, reltol = 1e-12)), error = function(e) NULL)
    if (is.null(o1)) next
    o2 <- tryCatch(optim(o1$par, nll, gr, method = "BFGS",
                         control = list(maxit = 300, reltol = 1e-14)), error = function(e) NULL)
    o  <- if (!is.null(o2) && o2$value < o1$value) o2 else o1
    if (is.null(best) || o$value < best$value) best <- o
  }
  c(p1, best$par)
}

## ---------------------------------------------------------------------------
## Data and the four test statistics
## ---------------------------------------------------------------------------

# Simulate per-stratum sufficient statistics (the "data"). The default truth is
# interest = 1, nuisances = 0, residual = 1; boundary_k_sims.R supplies
# psi_true with the nuisances just off zero.
bnd_sim <- function(A, N, psi_true = NULL) {
  if (is.null(psi_true)) psi_true <- c(1, rep(0, ncol(A) - 2), 1)
  sigT <- as.vector(A %*% psi_true)
  sigT * rchisq(length(N), df = N)
}

# WLD, LRT, SCR_C, SCR_U in one pass; mirrors all_stats() in coverage_engine.R.
# The full fit is seeded with the constrained null fit so l_full >= l_null and
# the LRT is nonnegative by construction. All four are evaluated at the true
# value p1 of the interest parameter, which is 1 in boundary_k_sims.R.
bnd_all_stats <- function(A, N, SS, k, p1 = 1) {
  lik <- bnd_mk_lik(A, N)
  p   <- k + 2
  lb_full <- c(rep(0, k + 1), 1e-8)
  lb_red  <- c(rep(0, k), 1e-8)
  mom_res <- max(1e-3, min(SS / N))

  # constrained null fit (interest fixed at 1); shared by LRT and SCR_C
  Ared <- A[, -1, drop = FALSE]; offr <- p1 * A[, 1]
  nll_red <- function(eta) { s <- offr + as.vector(Ared %*% eta)
    if (any(s <= 0)) 1e12 else 0.5 * sum(N * log(s) + SS / s) }
  gr_red  <- function(eta) { s <- offr + as.vector(Ared %*% eta)
    if (any(s <= 0)) rep(0, k + 1) else as.vector(crossprod(Ared, 0.5 * (N / s - SS / s^2))) }
  red_starts <- c(list(c(rep(0, k), mom_res), rep(1, k + 1)),
                  lapply(1:5, function(i) pmax(lb_red, rnorm(k + 1, 0.3, 0.4))))
  red     <- bnd_fit_min(nll_red, gr_red, lb_red, red_starts)
  psi_red <- c(p1, red$par)

  # full constrained MLE (seeded with the null fit)
  full_starts <- c(list(psi_red, c(p1, rep(0, k), 1), rep(1, p)),
                   lapply(1:6, function(i) pmax(lb_full, rnorm(p, 0.5, 0.5))))
  full <- bnd_fit_min(function(ps) lik$nll(ps, SS), function(ps) lik$gr(ps, SS),
                      lb_full, full_starts)

  lrt   <- max(0, 2 * (red$value - full$value))
  wld   <- (full$par[1] - p1)^2 / solve(lik$info(full$par))[1, 1]  # = Delta^2 * eff.info
  scr_c <- lik$gr(psi_red, SS)[1]^2 * solve(lik$info(psi_red))[1, 1]
  psi_u <- bnd_fit_unrestricted(A, N, SS, k, p1)
  scr_u <- lik$gr(psi_u, SS)[1]^2 * solve(lik$info(psi_u))[1, 1]

  c(WLD = wld, LRT = lrt, SCR_C = scr_c, SCR_U = scr_u)
}
