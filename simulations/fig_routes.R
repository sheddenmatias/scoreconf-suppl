library(tidyverse)
library(patchwork)
library(ggthemes)
library(here)

PDF <- TRUE

# RECONF_RESDIR and RECONF_FIGDIR redirect the input and output folders
# run_local.R sets them so a local run never overwrites the paper's files
res_dir <- Sys.getenv("RECONF_RESDIR", unset = here("results"))
out_dir <- Sys.getenv("RECONF_FIGDIR", unset = here("figures"))
source(here("R", "ggplot_theme_Publication.R"))

lt_scale  <- scale_linetype_stat()
clr_scale <- scale_colour_stat()

LEV  <- c("WLD", "LRT", "SCR_C", "SCR_U")
YLIM <- c(0.89, 1)   # the restricted score in (b) leaves the panel

###############################################################################
# Two routes to the boundary in one design, correlated random intercept and
# slope, testing the covariance: (a) the variances shrink with the correlation
# fixed; (b) the correlation goes to one with the variances fixed. Both panels
# use 20 groups of 10 observations. Shared vertical range.
###############################################################################

cov_dat <- readRDS(file.path(res_dir, "cov_dat.Rds")) %>%
  mutate(stat = factor(stat, levels = LEV))

rho_file <- file.path(res_dir, "covrho_dat.Rds")
if (!file.exists(rho_file))
  stop(rho_file, " is missing; run the covrho sweep first")
covrho <- readRDS(rho_file) %>%
  mutate(stat = factor(stat, levels = LEV))

panel <- function(dat, xvar, title, xlab, ylab, legend, xscale = NULL) {
  p <- ggplot(dat, aes(x = .data[[xvar]], y = prop, group = stat)) +
    geom_line(aes(color = stat, lty = stat)) +
    geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)),
                alpha = 0.1) +
    geom_hline(yintercept = 0.95) +
    coord_cartesian(ylim = YLIM) +
    labs(lty = "Statistic", color = "Statistic", x = xlab, y = ylab) +
    lt_scale + clr_scale + theme_Publication() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1))) +
    ggtitle(title)
  if (!is.null(xscale)) p <- p + xscale
  p
}

# Both panels on a natural-log axis with ticks in original units. (b) is
# plotted as the distance to singularity 1 - rho, so both parameters move
# away from the boundary as you move left to right.
covrho$x <- 1 - covrho$rho

pa <- panel(filter(cov_dat, type == "corr", n1 == 20, psi1 <= 1.5), "psi1",
            "(a) Variances to zero",
            expression(psi[1] == psi[3]), "Coverage", TRUE,
            xscale = scale_x_continuous(transform = "log",
                                        breaks = c(0.005, 0.02, 0.1, 0.5, 1.5),
                                        labels = c("0.005", "0.02", "0.1",
                                                   "0.5", "1.5")))

pb <- panel(covrho, "x",
            "(b) Correlation to one",
            expression(1 - rho), "", FALSE,
            xscale = scale_x_continuous(transform = "log",
                                        breaks = c(0.001, 0.01, 0.1, 1),
                                        labels = c("0.001", "0.01", "0.1", "1")))

p_all <- pa + pb + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

if (PDF) ggsave(filename = file.path(out_dir, "fig_cov_routes.pdf"),
                plot = p_all, width = 8, height = 4.8, create.dir = TRUE)

# Numbers for the caption.
covrho %>%
  group_by(stat) %>%
  summarize(min = min(prop), max = max(prop), .groups = "drop") %>%
  print()
