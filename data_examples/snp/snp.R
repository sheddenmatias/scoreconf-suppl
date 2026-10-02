## Partitioned heritability of white blood cell count in the Diversity Outbred
## mice: everything behind Table 2.  Produces snp.Rds.
##
## Autosomal SNPs are stratified into four bins by minor allele frequency,
## giving four random effects plus a residual.  Three methods are compared:
## profile score intervals, restricted profile likelihood intervals computed
## directly, and universal inference.
##
## Runtime about two hours.

suppressMessages({library(qtl2); library(reconf); library(lme4qtl); library(lme4)})

DATA <- "do.zip"
B_SPLITS <- 10L          # sample splits for universal inference
ALPHA <- 0.05

## ---- data and model ---------------------------------------------------------
do <- read_cross2(DATA)
dosage_raw <- do.call(cbind, lapply(do$geno[!do$is_x_chr], function(g){
  d <- matrix(NA_real_, nrow(g), ncol(g), dimnames = dimnames(g))
  d[g == 1] <- 0; d[g == 2] <- 1; d[g == 3] <- 2; d }))
maf <- apply(dosage_raw, 2, function(x){p <- mean(x, na.rm = TRUE)/2; min(p, 1 - p)})

build_grm <- function(D){
  for (j in seq_len(ncol(D))) {na <- is.na(D[, j]); if (any(na)) D[na, j] <- mean(D[, j], na.rm = TRUE)}
  p <- colMeans(D)/2; s <- sqrt(2*p*(1 - p))
  tcrossprod(sweep(sweep(D, 2, 2*p, "-"), 2, s, "/"))/ncol(D) }
norm_grm <- function(K) K*nrow(K)/sum(diag(K))

y <- do$pheno[, "WBC"]; names(y) <- rownames(do$pheno)
keep <- !is.na(y); y <- y[keep]; n <- length(y)
## Four minor-allele-frequency bins, following the stratification of Yang et al.
## (2015).
brk <- c(0, 0.05, 0.15, 0.3, 0.5)
sets <- split(names(maf), cut(maf, breaks = brk, right = FALSE, include.lowest = TRUE))
nmark <- sapply(sets, length)
cat("n =", n, "; markers per bin:", nmark, "\n")

## The relationship matrices are used as they come out of the markers, without
## a ridge; those from few markers are singular.
Ks <- lapply(sets, function(s){
  K <- norm_grm(build_grm(dosage_raw[keep, s, drop = FALSE]))
  rownames(K) <- colnames(K) <- names(y); K })

ids <- paste0("b", seq_along(Ks))
dat <- as.data.frame(setNames(lapply(ids, function(i) names(y)), ids), stringsAsFactors = FALSE)
dat$y <- as.numeric(y); dat$Sex <- do$covar[names(y), "Sex"]
dat$ngen <- do$covar[names(y), "ngen"]
bd <- as.Date(do$covar[names(y), "birth_date"]); dat$birth_days <- as.numeric(bd - min(bd))
fit <- relmatLmer(as.formula(paste("y ~ 1 + Sex + ngen + birth_days +",
       paste0("(1|", ids, ")", collapse = " + "))), data = dat,
       relmat = setNames(Ks, ids), verbose = 0L)
cat("lme4 isSingular:", isSingular(fit), "\n")

grm_cor <- cor(sapply(Ks, function(K) K[lower.tri(K)]))

## ---- likelihood pieces ------------------------------------------------------
Y <- getME(fit,"y"); X <- getME(fit,"X"); Z <- getME(fit,"Z")
Hl <- reconf:::get_Hlist_lmer(fit)
pc <- reconf:::get_precomp(Y, X, Z, REML = TRUE, Hlist = Hl)
psi_hat <- reconf:::get_psi_hat_lmer(fit)
val <- function(p) reconf:::loglikelihood(psi = p, Y=Y, X=X, Z=Z, Hlist=Hl,
        REML=TRUE, get_val=TRUE, get_score=FALSE, get_inf=FALSE,
        precomp=pc, check=FALSE)$value
scr <- function(p) reconf:::loglikelihood(psi = p, Y=Y, X=X, Z=Z, Hlist=Hl,
        REML=TRUE, get_val=FALSE, get_score=TRUE, get_inf=FALSE,
        precomp=pc, check=FALSE)$score
cat("psi_hat:", signif(psi_hat, 6), "\n")

## ---- kappa ------------------------------------------------------------------
## kappa_j is the multiple correlation between the score for component j and the
## scores for the rest; (1-kappa_j^2)^{-1} is the variance inflation the nuisance
## parameters cause for component j.
I <- reconf:::loglikelihood(psi = psi_hat, Y=Y, X=X, Z=Z, Hlist=Hl, REML=TRUE,
       get_val=FALSE, get_score=FALSE, get_inf=TRUE, precomp=pc, check=FALSE)$inf_mat
kappa <- sapply(seq_len(nrow(I)), function(j){
  a <- I[j, -j, drop=FALSE]
  sqrt(as.numeric(a %*% solve(I[-j, -j, drop=FALSE], t(a))) / I[j, j]) })
cat("kappa        :", sprintf("%7.4f", kappa), "\n")
cat("1/(1-kappa^2):", sprintf("%7.2f", 1/(1-kappa^2)), "\n")

## ---- profile score intervals ------------------------------------------------
t0 <- Sys.time()
score <- ci_all_lmer(fit, onestep = FALSE)
t_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("score interval time (s):", round(t_score, 1), "\n")

## ---- restricted profile likelihood intervals, computed directly -------------
## Fix a component, maximize over the rest subject to nonnegativity; L-BFGS-B
## with the analytic score, several starts.
LO <- 1e-10
prof <- function(rho, idx) {
  free <- setdiff(seq_along(psi_hat), idx)
  f  <- function(q){p <- numeric(5); p[idx] <- rho; p[free] <- q; -val(p)}
  gr <- function(q){p <- numeric(5); p[idx] <- rho; p[free] <- q; -scr(p)[free]}
  starts <- list(psi_hat[free], pmax(psi_hat[free], 0.05),
                 rep(mean(psi_hat), length(free)), psi_hat[free]*0.5 + 0.05)
  best <- -Inf
  for (s in starts) {
    o <- try(stats::optim(pmax(s, LO), f, gr, method = "L-BFGS-B",
                          lower = rep(LO, length(free)),
                          control = list(factr = 1e5, maxit = 500)), silent = TRUE)
    if (!inherits(o, "try-error") && -o$value > best) best <- -o$value
  }
  best
}
mx <- prof(psi_hat[1], 1)
stopifnot(abs(mx + REMLcrit(fit)/2) < 1e-4)   # our criterion == lme4's

CUT <- qchisq(1 - ALPHA, 1)
ends <- function(idx) {
  g <- function(rho) 2*(mx - prof(rho, idx)) - CUT
  hat <- psi_hat[idx]; step <- max(hat*0.3, 0.05)
  ## lower: either a root, or the boundary at zero when zero is inside
  lo <- if (g(LO) <= 0) 0 else {
    x <- hat; repeat {xn <- max(x - step, LO)
      if (g(xn) > 0) break; x <- xn; if (xn <= LO) break}
    if (g(max(x - step, LO)) > 0) stats::uniroot(g, c(max(x-step, LO), x), tol=1e-6)$root else 0 }
  x <- hat
  repeat {xn <- x + step; if (g(xn) > 0) break; x <- xn; if (x > 1e4) break}
  c(lo, stats::uniroot(g, c(x, x + step), tol = 1e-6)$root)
}
t0 <- Sys.time()
profile_ci <- t(vapply(seq_along(psi_hat), ends, numeric(2)))
t_prof <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("profile interval time (s):", round(t_prof, 1), "\n")

## Likelihood ratio at the endpoints.  Upper endpoints should equal CUT; a lower
## endpoint of zero is the constraint rather than a root, so its ratio is below.
lr_check <- t(vapply(seq_along(psi_hat), function(j)
  c(2*(mx - prof(profile_ci[j,1], j)), 2*(mx - prof(profile_ci[j,2], j))), numeric(2)))

## ---- universal inference ----------------------------------------------------
## Zhang, Ekvall and Molstad (2026), package varcomp,
## https://github.com/yqzhang5972/varcomp.  Skipped when varcomp is absent.
ui <- NULL
if (requireNamespace("varcomp", quietly = TRUE)) {
  CRIT <- 1/ALPHA
  ## The split likelihood ratio code assumes mean zero, so the fixed effects are
  ## removed by an orthonormal basis of the orthogonal complement of col(X);
  ## what remains is the restricted likelihood used elsewhere.
  qrX <- qr(X); p_rank <- qrX$rank
  Q <- qr.Q(qrX, complete = TRUE)[, (p_rank + 1):n, drop = FALSE]
  yt <- as.numeric(crossprod(Q, Y))
  Kt <- lapply(Ks, function(K) crossprod(Q, K %*% Q))
  N <- length(yt)
  ll1 <- varcomp:::loglik1_naive_s; ll01 <- varcomp:::loglik01_naive_s

  ## Estimate under the alternative from the first half alone.
  fit_alt <- function(i1) {
    y1 <- yt[i1]; v <- as.numeric(var(y1))
    starts <- list(rep(1, 5), rep(v/5, 5), c(rep(1e-6, 4), v),
                   c(rep(v/10, 4), v/2), c(rep(v/20, 4), 0.9*v))
    best <- NULL; bv <- -Inf
    for (s in starts) {
      o <- try(stats::optim(pmax(s, 1e-8), ll1, control = list(fnscale = -1),
                            y1 = y1, K1_list = Kt[1], K2_list = Kt[2:4], i1 = i1,
                            method = "L-BFGS-B", lower = rep(1e-8, 5)), silent = TRUE)
      if (!inherits(o, "try-error") && o$value > bv) {bv <- o$value; best <- o$par}
    }
    best
  }
  ## For components 1 to 4 this is varcomp::slrt_naive; for the residual the
  ## interest matrix is the identity, whose coefficient varcomp adds itself, so
  ## that coefficient is held at zero.
  ui_lr <- function(rho, idx, i1, i0, opt1) {
    y1 <- yt[i1]; y0 <- yt[i0]
    if (idx <= 4) {
      K1 <- Kt[idx]; K2 <- Kt[-idx]; nfree <- 5L
      starts <- list(rep(1, 5), psi_hat[-idx], c(rep(1e-6, 3), psi_hat[5]))
    } else {
      K1 <- list(diag(N)); K2 <- Kt; nfree <- 4L
      starts <- list(rep(1, 4), psi_hat[1:4], rep(1e-6, 4))
    }
    obj <- function(par) ll01(par = if (idx <= 4) par else c(par, 0), rho = rho,
                              y0 = y0, y1 = y1, K1_list = K1, K2_list = K2,
                              i1 = i1, i0 = i0)
    best <- -Inf
    for (s in starts) {
      o <- try(stats::optim(pmax(s[seq_len(nfree)], 1e-8), obj,
                            control = list(fnscale = -1), method = "L-BFGS-B",
                            lower = rep(1e-8, nfree)), silent = TRUE)
      if (!inherits(o, "try-error") && o$value > best) best <- o$value
    }
    as.numeric(ll01(par = opt1[-1], rho = opt1[1], y0 = y0, y1 = y1,
                    K1_list = Kt[1], K2_list = Kt[2:4], i1 = i1, i0 = i0) - best)
  }
  invert <- function(g, hi_start) {
    grid <- seq(0, hi_start, length.out = 25)
    gv <- vapply(grid, g, numeric(1)); jm <- which.min(gv)
    if (gv[jm] > 0) return(c(NA, NA))
    lo <- if (gv[1] <= 0) 0 else
      stats::uniroot(g, c(grid[max(which(gv > 0 & seq_along(gv) < jm))], grid[jm]), tol = 1e-4)$root
    if (any(gv > 0 & seq_along(gv) > jm)) {
      hi_r <- grid[min(which(gv > 0 & seq_along(gv) > jm))]
    } else {
      hi_r <- grid[length(grid)]
      while (g(hi_r) <= 0 && hi_r < 1e4) hi_r <- hi_r*2
      if (hi_r >= 1e4) return(c(lo, Inf))
    }
    c(lo, stats::uniroot(g, c(grid[jm], hi_r), tol = 1e-4)$root)
  }
  HI <- c(1.5, 3, 4, 5, 8)
  splits <- lapply(seq_len(B_SPLITS), function(sd){
    set.seed(sd); i1 <- sort(sample(N, floor(N/2)))
    list(i1 = i1, i0 = setdiff(seq_len(N), i1)) })
  ## Timed in three pieces: the alternative fits, the per-split array, and the
  ## averaged inversion.
  t0 <- Sys.time()
  opts <- lapply(splits, function(s) fit_alt(s$i1))
  t_ui_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0 <- Sys.time()
  per <- array(NA_real_, c(B_SPLITS, 5, 2))
  for (b in seq_len(B_SPLITS)) for (j in 1:5) {
    g <- function(rho) ui_lr(rho, j, splits[[b]]$i1, splits[[b]]$i0, opts[[b]]) - log(CRIT)
    per[b, j, ] <- invert(g, HI[j])
  }
  t_ui_per <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ## Averaged over the splits.
  t0 <- Sys.time()
  avg <- matrix(NA_real_, 5, 2, dimnames = list(c(names(sets), "Residual"),
                                                c("lower", "upper")))
  for (j in 1:5) {
    g <- function(rho) log(mean(exp(vapply(seq_len(B_SPLITS), function(b)
          ui_lr(rho, j, splits[[b]]$i1, splits[[b]]$i0, opts[[b]]), numeric(1))))) - log(CRIT)
    avg[j, ] <- invert(g, HI[j])
  }
  t_ui_avg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("universal time (s): fits", round(t_ui_fit, 1), "; averaged inversion",
      round(t_ui_avg, 1), "; per-split diagnostic", round(t_ui_per, 1), "\n")
  ui <- list(per = per, avg = avg, seeds = seq_len(B_SPLITS), N = N, p = p_rank,
             timing = c(fit = t_ui_fit, avg = t_ui_avg, per = t_ui_per))
} else {
  message("varcomp not installed: universal inference block skipped.")
}

## ---- save, then report ------------------------------------------------------
LAB <- c(names(sets), "Residual")
saveRDS(list(n = n, nmark = nmark, grm_cor = grm_cor, psi_hat = psi_hat,
             information = I, kappa = kappa, score = score,
             profile = profile_ci, profile_lr_check = lr_check, ui = ui,
             restricted_max = mx, labels = LAB,
             timing = c(score = t_score, profile = t_prof,
                        ## cost of the universal column as reported: the
                        ## alternative fits plus the averaged inversion.
                        universal = if (is.null(ui)) NA_real_
                                    else unname(ui$timing["fit"] + ui$timing["avg"]))),
        "snp.Rds")

sc <- as.data.frame(unclass(score))
cat("\n")
cat(sprintf("%-12s %10s %20s %20s %20s\n", "component", "estimate",
            "profile score", "profile likelihood", "universal"))
for (j in 1:5) {
  u <- if (is.null(ui)) c(NA, NA) else ui$avg[j, ]
  cat(sprintf("%-12s %10.3f (%8.3f,%8.3f) (%8.3f,%8.3f) (%8.3f,%8.3f)\n", LAB[j],
              psi_hat[j], sc$lower[j], sc$upper[j],
              profile_ci[j,1], profile_ci[j,2], u[1], u[2]))
}
cat("\nlikelihood ratio at profile endpoints (target", round(CUT, 4),
    "; a lower endpoint at zero is the constraint):\n")
for (j in 1:5) cat(sprintf("  %-12s %8.4f %8.4f\n", LAB[j], lr_check[j,1], lr_check[j,2]))
