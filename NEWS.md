# bootComb 0.1.0.9000 (current development version)

* Various bug fixes.
* Second simulation scenario added (for adjusted prevalence estimation situation).
* Ability to specify number of optimisation iterations for `identifyBetaPars`.
* `identifyBetaPars` is no longer an exported function; strictly internal
* added support for normal and Poisson distributions in addition the beta distribution
* added `getBetaFromCI`, a function that returns beta distribution functions from a given confidence interval (and in particular a random sampling function)
* added `getNormFromCI`, a function that returns normal distribution functions from a given confidence interval (and in particular a random sampling function)
* added `getPoisFromCI`, a function that returns Poisson distribution functions from a given confidence interval (and in particular a random sampling function)
* change the organisation of code within the source R code files

# bootComb 0.1.0

* First version to be released on CRAN.
