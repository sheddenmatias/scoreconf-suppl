## lme4 profile intervals for the autism socialization example: the lme4
## column of Table 1, on the variance and covariance scale.
## Produces autism_lme4.Rds.
##
## The fit is the weighted model of autism_weighted.R.  lme4 profiles the
## unrestricted likelihood even for a restricted fit, and where the back-spline
## of a profile fails, confint() reports the parameter's bound (0 or +-Inf) as
## the endpoint.  This script prints the failed splines and the stored bounds
## so those endpoints can be identified.
##
## The lme4 column of Table 1 in the paper is reproduced with lme4 2.0.1, 
## Matrix 1.7.5 and R 4.5.3.
##
## Run from this directory.

suppressMessages({library(lme4); library(WWGbook)})

data(autism, package = "WWGbook")
aut <- autism[!is.na(autism$vsae), ]
aut$sicd <- factor(aut$sicdegp, labels = c("Low", "Medium", "High"))
aut$wt   <- 1 / aut$age^2
f2 <- suppressMessages(suppressWarnings(
  lmer(vsae ~ age + sicd + (1 + age | childid), aut, REML = TRUE, weights = wt)))
grab <- function(expr) {
  w <- character()
  val <- withCallingHandlers(suppressMessages(expr),
          warning = function(x){w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")})
  list(value = val, warnings = unique(w))
}
## Each scale is profiled once and both intervals come off the profile object,
## so confint does not refit.
p_sd <- grab(profile(f2, which = "theta_", signames = FALSE))
p_vc <- grab(profile(f2, which = "theta_", signames = FALSE, prof.scale = "varcov"))
sdcor  <- list(ci = suppressWarnings(confint(p_sd$value)), warnings = p_sd$warnings)
varcov <- list(ci = suppressWarnings(confint(p_vc$value)), warnings = p_vc$warnings)
pv  <- p_vc$value
bad <- names(which(vapply(attr(pv, "backward"), inherits, NA, "error")))
cat("sd/correlation scale:\n"); print(sdcor$ci)
cat("variance/covariance scale:\n"); print(varcov$ci)
cat("back-splines that failed:", paste(bad, collapse = ", "), "\n")
cat("stored bounds, lower:", attr(pv, "lower"), " upper:", attr(pv, "upper"), "\n")
cat("warnings:\n"); for (m in varcov$warnings) cat("  -", m, "\n")

saveRDS(list(sdcor = sdcor, varcov = varcov, failed_splines = bad,
             bounds = list(lower = attr(pv, "lower"), upper = attr(pv, "upper")),
             session = sessionInfo()),
        "autism_lme4.Rds")
