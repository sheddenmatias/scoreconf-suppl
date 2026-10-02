## Autism socialization example, Section 5.1: Table 1 and
## Figure 4.  Produces autism_weighted.Rds.
##
## The model is fitted with known weights 1/age^2, taken as lmer prior weights.
## psi = (Var(int), Cov(int, slope), Var(slope), Residual).

suppressMessages({library(lme4); library(reconf); library(WWGbook)})

data(autism, package = "WWGbook")
aut <- autism[!is.na(autism$vsae), ]
aut$sicd <- factor(aut$sicdegp, labels = c("Low", "Medium", "High"))
aut$wt   <- 1 / aut$age^2

fit <- suppressWarnings(suppressMessages(
  lmer(vsae ~ age + sicd + (1 + age | childid), data = aut, REML = TRUE,
       weights = wt)))

n_obs <- nrow(aut); n_sub <- length(unique(aut$childid))

## ---- the observed heteroscedasticity, and what the weights do to it ---------
sd_raw <- tapply(aut$vsae, aut$age, sd)
sd_wtd <- tapply(residuals(fit) * sqrt(aut$wt), aut$age, sd)

## ---- likelihood pieces ------------------------------------------------------
## .lmer_matrices applies the prior weights: it returns W^(1/2) times each of
## Y, X and Z, which turns Var(E) = psi_r W^{-1} into Var(E) = psi_r I.
m  <- reconf:::.lmer_matrices(fit)
Hl <- reconf:::get_Hlist_lmer(fit)
pc <- reconf:::get_precomp(m$Y, m$X, m$Z, REML = TRUE, Hlist = Hl)
ll <- function(p) reconf:::loglikelihood(psi = p, Y = m$Y, X = m$X, Z = m$Z,
        Hlist = Hl, REML = TRUE, get_val = TRUE, get_score = FALSE,
        get_inf = FALSE, precomp = pc, check = FALSE)$value

psi_c <- reconf:::get_psi_hat_lmer(fit)
psi_u <- reconf:::maximize_loglik(start_val = psi_c, opt_idx = 1:4,
           Y = m$Y, X = m$X, Z = m$Z, Hlist = Hl, REML = TRUE, precomp = pc,
           check = FALSE, warn_nonconv = FALSE, iterlim = 5000L)$arg

## reconf's criterion equals lme4's up to sum(log(w))/2 from the weighting
stopifnot(abs(-2 * ll(psi_c) - REMLcrit(fit) - sum(log(aut$wt))) < 1e-6)

gap <- 2 * (ll(psi_u) - ll(psi_c))

Psi   <- matrix(c(psi_c[1], psi_c[2], psi_c[2], psi_c[3]), 2, 2)
eigPsi <- eigen(Psi)
rho   <- psi_c[2] / sqrt(psi_c[1] * psi_c[3])

## ---- kappa: multiple correlation between the score for component j and the
## scores for the rest.  (1 - kappa^2)^{-1} is the variance inflation the
## nuisance parameters cause for component j.
I <- reconf:::loglikelihood(psi = psi_c, Y = m$Y, X = m$X, Z = m$Z, Hlist = Hl,
       REML = TRUE, get_val = FALSE, get_score = FALSE, get_inf = TRUE,
       precomp = pc, check = FALSE)$inf_mat
kappa <- sapply(seq_len(nrow(I)), function(j) {
  a <- I[j, -j, drop = FALSE]
  sqrt(as.numeric(a %*% solve(I[-j, -j, drop = FALSE], t(a))) / I[j, j])
})

## variance each eigen-direction of Psi contributes to a weighted observation,
## relative to the error variance; the small one measures nearness to the cone
Zt      <- sqrt(aut$wt) * cbind(1, aut$age)
contrib <- sapply(1:2, function(k)
  eigPsi$values[k] * mean((Zt %*% eigPsi$vectors[, k])^2)) / psi_c[4]

## ---- profile score intervals ------------------------------------------------
## Section 5.1 reports the median of REPS runs for the two package timings.
REPS <- 5L
time_median <- function(f, reps = REPS) {
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    val <- f()
    ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  list(value = val, time = median(ts))
}

r_score <- time_median(function()
  ci_all_lmer(fit, test_idx = 1:4, level = 0.95, onestep = FALSE))
score   <- r_score$value
t_score <- r_score$time

## ---- restricted profile likelihood intervals, computed directly -------------
## Fix a component, maximize the restricted likelihood over the rest subject to
## positive semi-definiteness; feasibility is checked at every point.
psd <- function(p) p[1] >= 0 && p[3] >= 0 && p[4] > 0 && p[1]*p[3] - p[2]^2 >= -1e-10
prof_at <- function(val, idx) {
  free <- setdiff(1:4, idx)
  ## Nelder-Mead on the feasible set, four starts, each restarted from its own
  ## solution until it stops improving.
  obj <- function(q) {pp <- numeric(4); pp[idx] <- val; pp[free] <- q
                      if (!psd(pp)) -1e10 else ll(pp)}
  ## For a covariance far from its estimate the starting value must be inflated
  ## until Psi is positive semi-definite, or the simplex begins in the penalty.
  infl <- psi_c
  if (idx == 2 && val^2 > psi_c[1] * psi_c[3]) {
    s <- sqrt(val^2 / (psi_c[1] * psi_c[3])) * 1.05
    infl[1] <- psi_c[1] * s; infl[3] <- psi_c[3] * s
  }
  best <- -Inf
  for (s0 in list(infl[free], psi_c[free], psi_c[free] * 1.25, psi_c[free] * 0.8)) {
    par <- s0
    for (k in 1:10) {
      o <- try(stats::optim(par, obj, control = list(fnscale = -1, maxit = 5000,
             reltol = 1e-12), method = "Nelder-Mead"), silent = TRUE)
      if (inherits(o, "try-error")) break
      if (o$value <= best + 1e-10) {best <- max(best, o$value); break}
      best <- o$value; par <- o$par
    }
  }
  ## The unconstrained optimizer is faster and agrees when the constraint does
  ## not bind, which at this interior fit is the usual case; take the better.
  st <- psi_c; st[idx] <- val
  if (psd(st)) {
    o <- try(reconf:::maximize_loglik(start_val = st, opt_idx = free, Y = m$Y,
           X = m$X, Z = m$Z, Hlist = Hl, REML = TRUE, precomp = pc, check = FALSE,
           warn_nonconv = FALSE, iterlim = 5000L), silent = TRUE)
    if (!inherits(o, "try-error") && psd(o$arg)) best <- max(best, ll(o$arg))
  }
  best
}
mx  <- ll(psi_c)
CUT <- qchisq(0.95, 1)
lrt <- function(val, idx) 2 * (mx - prof_at(val, idx))
ends <- function(idx, lower_lim) {
  g <- function(v) lrt(v, idx) - CUT
  se <- sqrt(diag(solve(I)))[idx]
  x <- psi_c[idx]; step <- max(se, abs(x) * 0.1, 1e-6)
  lo <- NA; y <- x
  repeat {
    yn <- y - step
    if (!is.na(lower_lim) && yn <= lower_lim) {
      lo <- if (g(lower_lim + 1e-8) <= 0) lower_lim else
              stats::uniroot(g, c(lower_lim + 1e-8, y), tol = 1e-6)$root
      break }
    if (g(yn) > 0) {lo <- stats::uniroot(g, c(yn, y), tol = 1e-6)$root; break}
    y <- yn; if (y < -1e6) break }
  y <- x
  repeat {
    yn <- y + step
    if (g(yn) > 0) {hi <- stats::uniroot(g, c(y, yn), tol = 1e-6)$root; break}
    y <- yn; if (y > 1e6) {hi <- NA; break} }
  c(lo, hi)
}
LIM <- c(0, NA, 0, 0)
t0 <- Sys.time()
profile_ci <- t(vapply(1:4, function(j) ends(j, LIM[j]), numeric(2)))
t_prof <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
lr_check   <- t(vapply(1:4, function(j)
  c(lrt(profile_ci[j, 1], j), lrt(profile_ci[j, 2], j)), numeric(2)))

## ---- lme4, for the comparison column ---------------------------------------
## Default standard deviation and correlation scale.  parm = "theta_" asks for
## the four random-effect parameters only, so the timing covers the same number
## of intervals as the other methods.  The variance and covariance scale is in
## autism_lme4.R.
r_lme4 <- time_median(function()
  suppressWarnings(suppressMessages(
    confint(fit, method = "profile", oldNames = TRUE, parm = "theta_"))))
lme4_sdcor <- r_lme4$value
t_lme4     <- r_lme4$time

## ---- the three statistics as functions of the covariance, for the figure ----
## Wald uses the expected information at the fit.
se2  <- sqrt(solve(I)[2, 2])
grid <- seq(psi_c[2] - 3.2 * se2, psi_c[2] + 3.2 * se2, length.out = 61)
safe <- function(f) function(v) tryCatch(f(v), error = function(e) NA_real_)
## theta_null is the full parameter vector; test_idx says which element is held
## at its null value while the rest are re-estimated.
## theta_null also starts the nuisance maximization, so it is inflated the same
## way; the extended set only needs Sigma-tilde positive definite, but a start
## with Psi far from positive semi-definite can fail that.
score_at <- function(v) {
  tn <- psi_c; tn[2] <- v
  if (v^2 > psi_c[1] * psi_c[3]) {
    s <- sqrt(v^2 / (psi_c[1] * psi_c[3])) * 1.05
    tn[1] <- psi_c[1] * s; tn[3] <- psi_c[3] * s
  }
  as.numeric(unclass(suppressWarnings(
    score_test_lmer(fit, theta_null = tn, test_idx = 2)))[1])
}
curves <- data.frame(
  psi2  = grid,
  score = vapply(grid, safe(score_at), numeric(1)),
  lrt   = vapply(grid, safe(function(v) lrt(v, 2)), numeric(1)),
  wald  = ((grid - psi_c[2]) / se2)^2)
curves$lrt[curves$lrt > 1e6] <- NA_real_   # no feasible point at that value

## ---- save, then report ------------------------------------------------------
LAB <- c("Var(int)", "Cov(int,slope)", "Var(slope)", "Residual")
saveRDS(list(n_obs = n_obs, n_sub = n_sub, rho = rho, gap = gap,
             sd_raw = sd_raw, sd_weighted = sd_wtd,
             psi_constrained = psi_c, psi_unrestricted = psi_u,
             Psi = Psi, eigenvalues = eigPsi$values, contribution = contrib,
             information = I, kappa = kappa,
             score = score, profile = profile_ci, profile_lr_check = lr_check,
             lme4_sdcor = lme4_sdcor, curves = curves, cut = CUT,
             data = aut[, c("childid", "age", "vsae", "sicd", "wt")],
             labels = LAB),
        "autism_weighted.Rds")

## Everything Section 5.1 reports.  The lme4 column of Table 1
## is on the variance and covariance scale and comes from autism_lme4.R.
sc  <- as.data.frame(unclass(score))
fmt <- function(lo, hi) sprintf("(%.3f, %.3f)", lo, hi)
tab <- data.frame(component = LAB,
                  estimate  = round(psi_c, 3),
                  score     = fmt(sc$lower, sc$upper),
                  profile   = fmt(profile_ci[, 1], profile_ci[, 2]))
print(tab, row.names = FALSE)

cat(sprintf("\n%d observations on %d children\n", n_obs, n_sub))
cat(sprintf("response sd by age : %s  (factor %.1f)\n",
            paste(sprintf("%.1f", sd_raw), collapse = " "),
            max(sd_raw) / min(sd_raw)))
cat(sprintf("weighted resid sd  : factor %.1f\n", max(sd_wtd) / min(sd_wtd)))
cat(sprintf("correlation %.3f;  Psi eigenvalues %.1f %.2f;  contributions %.2f %.2f\n",
            rho, eigPsi$values[1], eigPsi$values[2], contrib[2], contrib[1]))
cat(sprintf("kappa              : %s\n", paste(sprintf("%.3f", kappa), collapse = " ")))
cat(sprintf("variance inflation : %s\n",
            paste(sprintf("%.0f", 1 / (1 - kappa^2)), collapse = " ")))
cat(sprintf("time (s), median of %d: score %.1f  lme4 %.1f\n",
            REPS, t_score, t_lme4))
cat("\nlme4 profile, standard deviation and correlation scale:\n")
print(lme4_sdcor)
