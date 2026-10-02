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

# Colors and line types for this figure, with the likelihood ratio limit
# (CONE) as a fifth series.
linetype_values <- c(SCR_U = "dashed", SCR_C = "dotted", LRT = "solid", WLD = "dotdash",
                     CONE = "longdash")
color_values    <- c(WLD = "#1a4e8c", LRT = "#b22222", SCR_C = "#1e6e1e", SCR_U = "#7b2d8b",
                     CONE = "grey45")
stat_labels["CONE"] <- "Likelihood ratio limit"

# r = k + 2: the interest parameter, k boundary nuisances, the error variance.
dat <- readRDS(file.path(res_dir, "boundary_k_dat_both.Rds")) %>% mutate(r = k + 2)

panel <- function(sn, title, ylab, legend) {
  d <- dat %>% filter(sense == sn)
  stat_dat <- d %>% filter(stat != "CONE") %>%
    mutate(stat = factor(stat, levels = c("WLD", "LRT", "SCR_C", "SCR_U")))
  cone_dat <- d %>% filter(stat == "CONE")
  ggplot(stat_dat, aes(x = r, y = prop, group = stat)) +
    geom_line(data = cone_dat, aes(color = stat, lty = stat)) +
    geom_point(data = cone_dat, aes(color = stat), shape = 4, size = 2) +
    geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
    geom_line(aes(color = stat, lty = stat)) +
    geom_hline(yintercept = 0.95) +
    scale_linetype_manual(values = linetype_values, labels = stat_labels,
                          breaks = c("WLD", "LRT", "SCR_C", "SCR_U", "CONE")) +
    scale_colour_manual(values = color_values, labels = stat_labels,
                        breaks = c("WLD", "LRT", "SCR_C", "SCR_U", "CONE")) +
    scale_x_continuous(breaks = sort(unique(stat_dat$r))) +
    coord_cartesian(ylim = range(dat$prop, na.rm = TRUE)) +
    labs(lty = "Statistic", color = "Statistic",
         x = "Number of variance parameters (r)", y = ylab) +
    theme_Publication() +
    theme(legend.position = if (legend) "bottom" else "none") +
    ggtitle(title)
}

pa <- panel("min", "(a) Least favorable", "Coverage", TRUE)
pb <- panel("max", "(b) Most favorable", "", FALSE)

fig <- pa + pb + plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")

if (PDF) ggsave(file.path(out_dir, "fig_cov_boundary_k_range.pdf"), fig,
                width = 9, height = 5, create.dir = TRUE)

dat %>%
  filter(stat %in% c("LRT", "CONE")) %>%
  select(sense, r, stat, prop) %>%
  pivot_wider(names_from = c(stat, sense), values_from = prop) %>%
  print(n = Inf)
