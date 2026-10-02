# Data generation for the simulations in Figures 1 and 2: the correlated
# random intercept and slope model, fitted by REML with lme4.

# Position of the tested component in lme4's ordering of the variance
# parameters. For the correlated random intercept and slope model lme4 uses
# the natural ordering (intercept variance, covariance, slope variance, error
# variance), so the position is idx itself.
grp_psi_idx <- function(fit, idx, type) {
  if (type == "corr") return(idx)
  stop("unknown simulation type: ", type)
}

# Correlated random intercept and slope: n1 groups of n2 observations.
sim_corr <- function(psi, n1, n2, p = 1) {
  psi0 <- psi[4]
  n    <- n1 * n2
  Psi1 <- matrix(c(psi[1], psi[2], psi[2], psi[3]), 2, 2)

  X_cov <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)
  colnames(X_cov) <- paste0("x", seq_len(p))
  X <- cbind(1, X_cov)

  # Z uses intercept + x1 for random slope
  Z <- Matrix::bdiag(lapply(seq_len(n1),
                            function(ii) X[((ii-1)*n2+1):(ii*n2), 1:2]))

  Psi1_eig <- eigen(Psi1)
  R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
  R  <- Matrix::kronecker(Matrix::Diagonal(n1), R1)

  beta  <- rnorm(p + 1)
  Xbeta <- drop(X %*% beta)
  y     <- Xbeta + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

  mc_data_frame <- data.frame(out = as.vector(y),
                               clust = as.factor(rep(seq_len(n1), each = n2)),
                               X_cov)

  cov_terms <- paste(paste0("x", seq_len(p)), collapse = " + ")
  fmla <- as.formula(paste("out ~", cov_terms, "+ (1 + x1|clust)"))
  lmer(fmla, data = mc_data_frame, REML = TRUE,
       control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
}

# Positions, within the psi vector, of a correlated two-by-two random-effect
# block (variance, covariance, variance). lme4 orders grouping factors by
# number of levels, so the block is not always first: a model with an extra
# independent effect on a finer factor puts that variance ahead of it. Returns
# NULL when every random effect is independent.
corr_block_idx <- function(fit) {
  vc <- lme4::VarCorr(fit)
  pos <- 0L
  for (nm in names(vc)) {
    q <- nrow(vc[[nm]])
    if (q == 2L) return(pos + c(1L, 2L, 3L))   # vech order: v11, v21, v22
    pos <- pos + q * (q + 1L) / 2L
  }
  NULL
}
