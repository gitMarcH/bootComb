## bootComb

This is an R package (under development) to provide functions to combine several estimates via a function to compute a new parameter.

Uncertainty is propagated through by sampling from probability distributions for each estimate and then deriving an empirical interval containing 95% of the resulting distribution (either via highest density interval or quantiles).

Several specific applications, though the methodology is quite general:
1. Combining several conditional prevalences (see e.g. [Stockdale et al., J. Hepatol. (2020)](https://doi.org/10.1016/j.jhep.2020.04.008))
2. Adjusting a prevalence for assay sensitivity and specificity when only estimates of the latter are known.
