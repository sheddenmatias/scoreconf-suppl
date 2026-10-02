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

# The one-step interval appears alongside the four shared statistics, extend
# the shared scales with a fifth color.
stat_colours["SCR_1"]   <- "#56B4E9"
stat_linetypes["SCR_1"] <- "longdash"
stat_labels["SCR_1"]    <- "Score, one-step"

###############################################################################
# Interval widths, from width_dat.Rds. The figure shows the width
# distributions of the likelihood-ratio and proposed intervals in the
# high-signal cells (snr_*), with the fraction of unbounded intervals, if any,
# printed above each violin. The script also prints the median paired width
# ratios quoted in Section 4 and, from the rho_0.999 cell, the fraction of
# empty constrained-score intervals.
###############################################################################

# The cells reported in the paper (see code/width_engine.R). The saved
# width_dat.Rds may hold more cells; only these are plotted and printed.
PAPER_CELLS <- c("snr_20", "snr_40", "snr_80", "snr_160", "rho_0.999")

dat <- readRDS(file.path(res_dir, "width_dat.Rds")) %>%
  filter(key %in% PAPER_CELLS) %>%
  mutate(width = upper - lower,
         unbounded = is.infinite(width))

lrt_w <- dat %>%
  filter(stat == "LRT") %>%
  select(key, rep, block, w_lrt = width)

ratios <- dat %>%
  filter(stat %in% c("SCR_U", "SCR_1")) %>%
  inner_join(lrt_w, by = c("key", "rep", "block")) %>%
  group_by(key, stat) %>%
  summarize(ratio = median(width / w_lrt, na.rm = TRUE), .groups = "drop")

## --- figure: width distributions in the high-signal cells -------------

dists <- dat %>%
  filter(grepl("^snr_", key), stat %in% c("LRT", "SCR_U")) %>%
  mutate(n1 = as.integer(sub("snr_", "", key)))

unb <- dists %>%
  group_by(n1, stat) %>%
  summarize(frac = mean(unbounded), top = max(width[is.finite(width)]),
            .groups = "drop") %>%
  filter(frac > 0)

p_dist <- ggplot(filter(dists, is.finite(width)),
                 aes(x = factor(n1), y = width, fill = stat)) +
  geom_violin(position = position_dodge(width = 0.8), linewidth = 0.2) +
  geom_text(data = unb, position = position_dodge2(width = 0.8),
            aes(x = factor(n1), y = top * 1.5, group = stat,
                label = sprintf("%.*f%% unb.", ifelse(frac < 0.01, 1, 0),
                                100 * frac)),
            size = 2.8, inherit.aes = FALSE) +
  scale_y_continuous(transform = "log", breaks = c(0.5, 1, 2, 4, 8)) +
  scale_fill_manual(values = stat_colours, labels = stat_labels,
                    name = "Statistic") +
  labs(x = "Number of groups", y = "Interval width") +
  theme_Publication() +
  theme(legend.position = "bottom")

if (PDF) ggsave(file.path(out_dir, "fig_width_dists.pdf"), p_dist,
                width = 8, height = 4.8, create.dir = TRUE)

## --- numbers for the captions and the text ---------------------------------

cover <- dat %>%
  mutate(cover = !is.na(lower) & lower <= true & true <= upper) %>%
  group_by(key, stat) %>%
  summarize(prop = mean(cover), unb = mean(unbounded), .groups = "drop")

cat("median width ratios to the likelihood ratio:\n")
ratios %>% pivot_wider(names_from = stat, values_from = ratio) %>%
  arrange(key) %>% print(n = Inf)
cat("\ncoverage and unbounded fraction:\n")
cover %>% filter(stat != "WLD") %>%
  pivot_wider(names_from = stat, values_from = c(prop, unb)) %>%
  arrange(key) %>% print(n = Inf)

cat("\nempty constrained-score (SCR_C) intervals at correlation 0.999:\n")
dat %>%
  filter(key == "rho_0.999", stat == "SCR_C") %>%
  summarize(empty = mean(!is.na(width) & width <= 0), n = n()) %>%
  print()
