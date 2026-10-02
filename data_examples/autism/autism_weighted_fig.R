## Figure 4 for the autism socialization example.  Reads
## autism_weighted.Rds and writes fig_autism.pdf next to it.
##
## Panel (a): observed trajectories with the mean and one standard deviation at
## each age.  Panel (b): each statistic as a function of the covariance, with
## the chi-squared cutoff and the interval each one returns.

suppressMessages({library(ggplot2); library(patchwork)})

a  <- readRDS("autism_weighted.Rds")
d  <- a$data
cu <- a$curves

## Okabe-Ito palette, with line type repeating the identity.
COL <- c(score = "#CC79A7", lrt = "#D55E00", wald = "#0072B2")
LTY <- c(score = "dashed",  lrt = "solid",   wald = "dotdash")
LAB <- c(score = "Score, extended", lrt = "Likelihood ratio", wald = "Wald")

suppressMessages(source("../../R/ggplot_theme_Publication.R"))
base_thm <- theme_Publication(base_size = 12)

## ggplot2 moved the inside-panel legend to its own argument at 3.5.0.
leg_inside <- function(pos) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = pos)
  else theme(legend.position = pos)
}

## ---- panel (a): trajectories, with the mean and one standard deviation ------
sumry <- data.frame(
  age  = as.numeric(names(tapply(d$vsae, d$age, mean))),
  mean = as.numeric(tapply(d$vsae, d$age, mean)),
  sd   = as.numeric(tapply(d$vsae, d$age, sd)))

pa <- ggplot() +
  geom_line(data = d, aes(x = age, y = vsae, group = childid),
            colour = "grey65", linewidth = 0.25, alpha = 0.6) +
  geom_errorbar(data = sumry, aes(x = age, ymin = mean - sd, ymax = mean + sd),
                width = 0.45, linewidth = 0.7, colour = "black") +
  geom_point(data = sumry, aes(x = age, y = mean), size = 2, colour = "black") +
  scale_x_continuous(breaks = c(2, 3, 5, 9, 13)) +
  labs(x = "Age (years)", y = "Socialization score",
       title = "(a) Observed trajectories") +
  base_thm

## ---- panel (b): the three statistics, and the intervals they give -----------
long <- rbind(
  data.frame(psi2 = cu$psi2, stat = "score", value = cu$score),
  data.frame(psi2 = cu$psi2, stat = "lrt",   value = cu$lrt),
  data.frame(psi2 = cu$psi2, stat = "wald",  value = cu$wald))
long <- long[is.finite(long$value), ]
long$stat <- factor(long$stat, levels = c("score", "lrt", "wald"))

## endpoints: where each curve crosses the cutoff
cross <- function(y) {
  i <- which(diff(sign(y - a$cut)) != 0)
  vapply(i, function(k) stats::approx(y[k:(k+1)], cu$psi2[k:(k+1)],
                                      xout = a$cut)$y, numeric(1))
}
ends <- rbind(
  data.frame(stat = "score", lo = a$score[2, 2], hi = a$score[2, 3]),
  data.frame(stat = "lrt",   lo = cross(cu$lrt)[1],  hi = cross(cu$lrt)[2]),
  data.frame(stat = "wald",  lo = cross(cu$wald)[1], hi = cross(cu$wald)[2]))
ends$stat <- factor(ends$stat, levels = c("score", "lrt", "wald"))
ytop <- max(long$value)
ends$y <- ytop * c(-0.10, -0.17, -0.24)

## The interval bars repeat their curve's line type.
pb <- ggplot() +
  geom_hline(yintercept = a$cut, colour = "grey40", linewidth = 0.4) +
  annotate("text", x = max(cu$psi2), y = a$cut, vjust = 1.6, hjust = 1,
           label = "chi[1]^2~0.95~quantile", parse = TRUE,
           size = 3.1, colour = "grey30") +
  geom_line(data = long, aes(x = psi2, y = value, colour = stat, linetype = stat),
            linewidth = 0.8) +
  geom_segment(data = ends, aes(x = lo, xend = hi, y = y, yend = y,
               colour = stat, linetype = stat), linewidth = 0.9,
               lineend = "butt", show.legend = FALSE) +
  geom_point(data = data.frame(x = a$psi_constrained[2], y = 0), aes(x = x, y = y),
             size = 2, colour = "black") +
  scale_colour_manual(values = COL, labels = LAB, name = NULL) +
  scale_linetype_manual(values = LTY, labels = LAB, name = NULL) +
  labs(x = "Covariance of intercept and slope", y = "Test statistic",
       title = "(b) Inverting the test") +
  coord_cartesian(ylim = c(ytop * -0.30, ytop)) +
  base_thm + leg_inside(c(0.30, 0.90)) +
  theme(legend.background = element_rect(fill = "white", colour = NA),
        legend.key.width = unit(1.5, "cm"))

p <- pa + pb + plot_layout(ncol = 2)
ggsave("fig_autism.pdf", p, width = 10, height = 4.6)
