# Pool the blocks written by run_cell.R into the files the figure scripts read.
#
# A block reports, for each statistic, the proportion covered and the number of
# replications that produced a value. Pooling is the weighted mean over blocks
# with the replication counts as weights, which is exactly the proportion the
# single long run would have given; the standard error is recomputed from the
# pooled proportion and the total count.
#
# Reports any block that is missing before writing anything, so a partial array
# is visible rather than silently producing a figure with thin cells.

root <- normalizePath(Sys.getenv("RECONF_DIR", unset = getwd()))
setwd(root)

# Flat on a cluster, simulations/code/ in the repository; see run_cell.R.
locate <- function(name) {
  for (p in c(file.path(root, name),
              file.path(root, "R", name),
              file.path(root, "simulations", "code", name)))
    if (file.exists(p)) return(p)
  stop("cannot find ", name, " under ", root)
}
here <- function(...) locate(basename(file.path(...)))

source(locate("cells.R"))

# RECONF_CELLDIR and RECONF_RESDIR match the overrides in run_cell.R, so a
# smoke test cannot overwrite the results the figures are built from.
cell_dir <- Sys.getenv("RECONF_CELLDIR",
                       unset = file.path(root, "results", "cells"))
res_dir  <- Sys.getenv("RECONF_RESDIR", unset = file.path(root, "results"))
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(cell_dir, pattern = "\\.Rds$", full.names = TRUE)
if (!length(files)) stop("no block files in ", cell_dir)

blocks <- lapply(files, readRDS)
all_blocks <- dplyr::bind_rows(blocks)
all_blocks <- all_blocks[all_blocks$key %in% names(CELLS), , drop = FALSE]
if (!nrow(all_blocks)) stop("no blocks in ", cell_dir, " belong to cells.R")

## --- completeness ---------------------------------------------------------
# Only sweeps with blocks on disk are checked. A deliberate single-sweep run
# (RECONF_SWEEP) should not be reported as missing the other one.
# RECONF_NSIM and RECONF_BLOCK must match the values the blocks were run with.
env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(v) || v < 1L) default else v
}
block_sz <- env_int("RECONF_BLOCK", 250L)
have    <- table(all_blocks$key[!duplicated(all_blocks[c("key", "block")])])
in_play <- CELLS[vapply(CELLS, `[[`, character(1), "sweep") %in%
                   unique(all_blocks$sweep)]
want <- vapply(in_play, function(x)
  ceiling(env_int("RECONF_NSIM", x$n_sim) / block_sz), numeric(1))
short <- names(want)[!names(want) %in% names(have) |
                     have[names(want)] < want]
short <- short[!is.na(short)]
if (length(short)) {
  cat("INCOMPLETE cells (", length(short), " of ", length(want), "):\n", sep = "")
  for (k in short)
    cat(sprintf("  %-28s %s of %s blocks\n", k,
                if (k %in% names(have)) have[[k]] else 0, want[[k]]))
} else {
  cat("all", length(want), "cells complete\n")
}

## --- pool -----------------------------------------------------------------
id_cols <- setdiff(names(all_blocks), c("prop", "se", "reps", "block"))

pooled <- all_blocks |>
  dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) |>
  dplyr::summarise(prop = sum(prop * reps) / sum(reps),
                   reps = sum(reps), .groups = "drop") |>
  dplyr::mutate(se = sqrt(prop * (1 - prop) / reps)) |>
  as.data.frame()

write_sweep <- function(sweep, file, cols = NULL) {
  d <- pooled[pooled$sweep == sweep, , drop = FALSE]
  if (!nrow(d)) {cat("skipping", file, "- no cells\n"); return(invisible(NULL))}
  d <- d[, setdiff(names(d), c("sweep", "key")), drop = FALSE]
  d <- d[, colSums(!is.na(d)) > 0, drop = FALSE]
  if (!is.null(cols)) d <- d[, cols, drop = FALSE]
  saveRDS(d, file.path(res_dir, file))
  cat("wrote", file, "-", nrow(d), "rows\n")
}

write_sweep("cov",    "cov_dat.Rds")
write_sweep("covrho", "covrho_dat.Rds")
