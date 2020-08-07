# bootComb

## Description

This is an R package to combine several estimates via a arbitrary function to compute a new parameter.

Uncertainty is propagated through by sampling from probability distributions for each estimate and then deriving an empirical interval containing 95% of the resulting distribution (either via highest density interval or quantiles).

Several specific applications, though the methodology is quite general:
1. Combining several conditional prevalences (see e.g. [Stockdale et al., J. Hepatol. (2020)](https://doi.org/10.1016/j.jhep.2020.04.008))
2. Adjusting a prevalence for assay sensitivity and specificity when the latter are not known exactly and only estimates are available.

The method assumes that the combined estimates are independent. This is generally the case where parameters estimated on different datasets are combined, but in some situations this may not hold. Future version of this package will include support for some joint distributions.

## Installation

### From CRAN (stable version)

\code{
install.packages('bootComb')
}

### From GitHub (development version)

\code{
# install.packages("devtools")
devtools::install_github("gitMarcH/bootComb")
}

## Example

Product of 2 probability parameters with beta distributions.

\code{
library(bootComb)

dist1<-function(n){rbeta(n,shape1=2,shape2=2)}
dist2<-function(n){rbeta(n,shape1=8,shape2=4)}

distListEx<-list(dist1,dist2)
combFunEx<-function(pars){pars[[1]]*pars[[2]]}

bootComb(distList=distListEx,combFun=combFunEx)
#> 
#> $conf.int
#>      lower      upper 
#> 0.03643913 0.64242678 
#> attr(,"credMass")
#> [1] 0.95
#> 
#> $bootstrapValues
#> NULL
}
