# Supplementary material

Code and data to reproduce results in "Score-based confidence intervals for 
variance-covariance parameters in linear mixed models" by Matias Shedden
and Karl Oskar Ekvall.

The accompanying R package `reconf` is available at
<https://github.com/koekvall/reconf>.

## Directories

- `R/` - Plot theme shared by the figure scripts
- `data_examples/snp/` - Mouse SNP example
- `data_examples/autism/` - Autism socialization example
- `simulations/` - Simulations, see `simulations/README.md`
- `results/` - generated `.Rds` from the simulations
- `figures/` - generated `.pdf` figures from the simulations

## Requirements

R packages:

- **Simulations** (`simulations/`): `reconf`, `Matrix`, `lme4`, `alabama`,
  `nloptr`, `foreach`, `doParallel`, `doRNG`, `dplyr`, `here`. The
  simulations used `reconf` commit `164d7bd`.
- **Autism example** (`data_examples/autism/`): `reconf`, `lme4`, `WWGbook`.
- **SNP example** (`data_examples/snp/`): `reconf`, `lme4`, `qtl2`, `lme4qtl`,
  and `varcomp` (<https://github.com/yqzhang5972/varcomp>) for the universal
  inference column of Table 2.
- **Figures**: `tidyverse` with `ggplot2` 3.5.0 or later, `patchwork`,
  `ggthemes`, `here`.
