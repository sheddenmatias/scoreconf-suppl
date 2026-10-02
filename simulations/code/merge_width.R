# Pool the width blocks into width_dat.Rds: one row per replication and
# statistic, with interval endpoints. Reports any block that is missing
# before writing anything, and the rate of failed intervals by cell and
# statistic.

root <- normalizePath(Sys.getenv("RECONF_DIR", unset = getwd()))
setwd(root)

locate <- function(name) {
  for (p in c(file.path(root, name),
              file.path(root, "R", name),
              file.path(root, "simulations", "code", name)))
    if (file.exists(p)) return(p)
  stop("cannot find ", name, " under ", root)
}
here <- function(...) locate(basename(file.path(...)))

source(locate("width_engine.R"))

env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name)))
  if (is.na(v) || v < 1L) default else v
}
block_sz <- env_int("RECONF_BLOCK", 250L)
n_sim    <- env_int("RECONF_NSIM", 1000L)

width_dir <- Sys.getenv("RECONF_WIDTHDIR",
                        unset = file.path(root, "results", "width_blocks"))
res_dir   <- Sys.getenv("RECONF_RESDIR", unset = file.path(root, "results"))
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(width_dir, pattern = "\\.Rds$", full.names = TRUE)
if (!length(files)) stop("no block files in ", width_dir)

dat <- do.call(rbind, lapply(files, readRDS))

cells <- width_cells(n_sim)
want  <- ceiling(vapply(cells, `[[`, numeric(1), "n_sim") / block_sz)
have  <- table(dat$key[!duplicated(dat[c("key", "block")])])
short <- names(want)[!names(want) %in% names(have) |
                     have[names(want)] < want]
short <- short[!is.na(short)]
if (length(short)) {
  cat("INCOMPLETE cells (", length(short), " of ", length(want), "):\n", sep = "")
  for (k in short)
    cat(sprintf("  %-14s %s of %s blocks\n", k,
                if (k %in% names(have)) have[[k]] else 0, want[[k]]))
} else {
  cat("all", length(want), "cells complete\n")
}

saveRDS(dat, file.path(res_dir, "width_dat.Rds"))
cat("wrote width_dat.Rds:", nrow(dat), "rows\n")

cat("\n=== interval failure rates by cell and statistic ===\n")
dat$failed <- is.na(dat$lower) | is.na(dat$upper)
print(aggregate(failed ~ key + stat, dat, mean), row.names = FALSE)
