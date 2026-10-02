# Autism socialization

Section 5.1, Table 1 and Figure 4.

## Reproduce

Run from this directory, in order:

```r
source("autism_weighted.R")      # writes autism_weighted.Rds
source("autism_weighted_fig.R")  # reads it, writes fig_autism.pdf
source("autism_lme4.R")          # lme4 column of the table; writes autism_lme4.Rds
```

`autism_weighted.R` prints the estimates and the profile-score and
profile-likelihood intervals. It also times `lme4`'s profile intervals on
lme4's default standard-deviation and correlation scale. It takes about two
minutes, most of it the five timing repetitions whose median Section 5.1
reports. The figure script takes a few seconds.

`autism_lme4.R` computes the `lme4` column of the table. It prints lme4's
profile intervals on the variance and covariance scale, the profiles whose
back-spline failed, and the parameter bounds that `confint` reports in their
place. It does not depend on the other two scripts.

## Data

`WWGbook::autism`, from CRAN. 158 children, scores at ages 2, 3, 5, 9 and 13,
610 observations after dropping missing responses.

## Requirements

`reconf`, `lme4`, `WWGbook`, `ggplot2`, `patchwork`. The figure script reads
the plot theme from `../../R/ggplot_theme_Publication.R`.
