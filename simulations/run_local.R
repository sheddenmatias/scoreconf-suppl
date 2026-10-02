## Run the simulations behind one figure of the paper, with the paper's settings
## and seeds, and draw that figure from the new results.
##
## Set the four variables below, then source this file. 
## Results go to local_runs/<figure>_nsim<N>/ in the repository root.

## ---- options -----------------------------------------------------------------

FIGURE <- "routes"   # "routes"     Figure 1 (coverage, two routes to the boundary)
                     # "width"      Figure 2 (interval widths)
                     # "boundary_k" Figure 3 (number of variance parameters)
                     # "all"        all three
NSIM   <- 200        # replications per setting. NA for the paper's (10000 for
                     # routes and boundary_k, 1000 for width), which takes long
CORES  <- NA         # parallel cores. NA for all cores but one
FRESH  <- FALSE      # TRUE deletes this run's output folder first

## -----------------------------------------------------------------------------

PAPER_NSIM <- c(routes = 10000L, width = 1000L, boundary_k = 10000L)
if (!FIGURE %in% c(names(PAPER_NSIM), "all"))
  stop('FIGURE must be "routes", "width", "boundary_k" or "all"')
if (!is.na(NSIM)  && NSIM  < 1) stop("NSIM must be at least 1, or NA")
if (!is.na(CORES) && CORES < 1) stop("CORES must be at least 1, or NA")

## The repository root is the parent of the folder holding this script: taken
## from source() when the file is sourced, from the RStudio editor when it is
## run line by line, and otherwise from the working directory.
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(this_file) && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable())
  this_file <- rstudioapi::getSourceEditorContext()$path
root <- if (length(this_file) && nzchar(this_file)) {
  normalizePath(file.path(dirname(this_file), ".."))
} else getwd()
if (!file.exists(file.path(root, "simulations", "code", "cells.R")))
  stop("cannot find the repository root; setwd() to it and source this file again")

## Each step runs in a fresh R process started from this R installation, with
## this session's package libraries, so it uses the same R and packages as the
## session that sources this file.
rscript <- file.path(R.home("bin"), "Rscript")

## The scripts take the repository root from the working directory; set it for
## each step and restore it afterwards.
in_dir <- function(dir, expr) {
  old <- setwd(dir)
  on.exit(setwd(old))
  expr
}

run <- function(script, env) {
  env <- c(env, R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  old <- Sys.getenv(names(env), unset = NA)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    for (nm in names(env))
      if (is.na(old[[nm]])) Sys.unsetenv(nm)
      else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
  })
  cat("\n==>", script, "\n")
  status <- in_dir(root, system2(rscript, shQuote(file.path(root, script))))
  if (!identical(status, 0L)) stop(script, " failed (exit status ", status, ")")
}

## Each figure: the simulation scripts, then the figure script.
steps <- list(
  routes = list(
    env   = function(dir) c(RECONF_SWEEP   = "cov,covrho",
                            RECONF_CELLDIR = file.path(dir, "blocks")),
    sims  = c("simulations/code/run_cell.R", "simulations/code/merge_cells.R"),
    fig   = "simulations/fig_routes.R",
    pdf   = "fig_cov_routes.pdf"),
  width = list(
    env   = function(dir) c(RECONF_WIDTHDIR = file.path(dir, "blocks")),
    sims  = c("simulations/code/run_width.R", "simulations/code/merge_width.R"),
    fig   = "simulations/fig_width.R",
    pdf   = "fig_width_dists.pdf"),
  boundary_k = list(
    env   = function(dir) character(0),
    sims  = "simulations/code/boundary_k_sims.R",
    fig   = "simulations/fig_boundary_k.R",
    pdf   = "fig_cov_boundary_k_range.pdf")
)

run_figure <- function(fig) {
  n   <- if (is.na(NSIM)) PAPER_NSIM[[fig]] else as.integer(NSIM)
  dir <- file.path(root, "local_runs", sprintf("%s_nsim%d", fig, n))
  if (FRESH && dir.exists(dir)) unlink(dir, recursive = TRUE)
  dir.create(file.path(dir, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir)

  env <- c(RECONF_NSIM   = as.character(n),
           RECONF_RESDIR = file.path(dir, "results"),
           RECONF_FIGDIR = file.path(dir, "figures"),
           steps[[fig]]$env(dir))
  if (!is.na(CORES)) env <- c(env, RECONF_CORES = as.character(as.integer(CORES)))

  cat(sprintf("\n######## %s: %d replications per setting, output in %s\n",
              fig, n, dir))
  for (s in steps[[fig]]$sims) run(s, env)
  run(steps[[fig]]$fig, env)
  cat(sprintf("\n######## %s done: %s\n", fig,
              file.path(dir, "figures", steps[[fig]]$pdf)))
}

for (fig in if (FIGURE == "all") names(PAPER_NSIM) else FIGURE) run_figure(fig)
