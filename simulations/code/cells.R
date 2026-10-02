# The 20 simulation settings behind Figure 1. Each setting (cell) is one
# value on a horizontal axis of the figure: a true parameter psi, simulated
# 10,000 times, giving the coverage of the nominal 95% interval for each of the
# four statistics (Wald, likelihood ratio, constrained score, proposed score).
# The settings are listed in a fixed order, so that a SLURM array index
# identifies one. Sourced by run_cell.R and merge_cells.R.
#
# Both sweeps use the correlated random intercept and slope model with 20
# groups of 10 observations and error variance psi4 = 1. The parameter of
# interest is the covariance psi2 (idx = 2).
#
#   cov     panel (a): psi1 = psi3 shrink toward zero, correlation fixed at 0.2
#   covrho  panel (b): psi1 = psi3 = 2.5, correlation rho goes to one
#
# Each cell has its own seed, and run_cell.R seeds block b of a cell with
# seed + b. A cell's results therefore do not depend on its position in this
# table or on which other cells are run. The keys and seeds are those of the
# original HiPerGator run.

N_SIM <- 10000L

# One entry per cell. `args` goes to run_coverage_sim; `extra` holds columns the
# figure scripts expect but the engine does not return; `seed` is the base for
# this cell, offset by the block index inside run_cell.R.
cell_list <- list()
add_cell <- function(sweep, key, args, extra = list(), seed) {
  cell_list[[length(cell_list) + 1L]] <<-
    list(sweep = sweep, key = key, args = args, extra = extra,
         seed = seed, n_sim = N_SIM)
}

## --- cov: panel (a), the variances shrink ----------------------------------
# Ten values of psi1 = psi3, roughly log-spaced from 0.005 to 1.5. The seeds
# were assigned in two batches in the original run, hence the two seed bases.
corr_fixed <- 0.2
cov_psi  <- c(0.005,   0.01,    0.02,    0.05,    0.1,
              0.15,    0.25,    0.5,     1,       1.5)
cov_seed <- c(102001L, 109102L, 102002L, 109101L, 102003L,
              109103L, 102004L, 102005L, 102006L, 109104L)

for (i in seq_along(cov_psi)) {
  pd <- cov_psi[i]
  add_cell("cov", sprintf("cov_corr_small_%g", pd),
           list(psi = c(pd, corr_fixed * pd, pd, 1),
                n1 = 20, n2 = 10, type = "corr", idx = 2),
           seed = cov_seed[i])
}

## --- covrho: panel (b), the correlation goes to one ------------------------
# Ten values of rho, with 1 - rho log-spaced from 1 down to 0.001.
psi_diag_fixed <- 2.5
rho_grid <- c(0,       0.5,     0.8,     0.9,     0.95,
              0.98,    0.99,    0.995,   0.998,   0.999)
rho_seed <- c(106001L, 106002L, 106003L, 106004L, 106005L,
              109201L, 106006L, 109202L, 109203L, 106007L)

for (i in seq_along(rho_grid)) {
  rho_i <- rho_grid[i]
  add_cell("covrho", sprintf("covrho_%g", rho_i),
           list(psi = c(psi_diag_fixed, rho_i * psi_diag_fixed,
                        psi_diag_fixed, 1),
                n1 = 20, n2 = 10, type = "corr", idx = 2),
           extra = list(rho = rho_i),
           seed = rho_seed[i])
}

CELLS <- cell_list
names(CELLS) <- vapply(CELLS, `[[`, character(1), "key")

if (anyDuplicated(names(CELLS))) stop("duplicate cell keys")
