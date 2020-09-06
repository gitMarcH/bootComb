# bootComb 0.2.0.9000 (current development version)

* Added cranlogs badges to README.md
* Added support fro negative binomial and gamma distributions in addition to the beta, normal and Poisson distributions.

?? remove suppot for Poisson and Neg Binomial? Quantile matching does not work as optim requires continuous functions to optimise, yet here the quantile functions are step functions...

# bootComb 0.2.0 (first CRAN version)

* Various bug fixes.
* Second simulation scenario added (for adjusted prevalence estimation situation).
* Ability to specify number of optimisation iterations for `identifyBetaPars`
* `identifyBetaPars` is no longer an exported function; strictly internal
* Added `getBetaFromCI`, a function that returns beta distribution functions from a given confidence interval (and in particular a random sampling function)
* Added support for normal and Poisson distributions in addition to the beta distribution
* Changed the organisation of code within the source R code files
* Implemented changes requested from the submission of v0.1.0 to CRAN (references in DESCRIPTION file, TRUE/FALSE instead of T/F, restoring user parameters after doing plots)

# bootComb 0.1.0

* First version to be submitted to CRAN (but rejected).
