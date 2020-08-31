# bootComb 0.2.0.9000 (current development version)

* Added cranlogs badges to README.md

# bootComb 0.2.0

* Various bug fixes.
* Second simulation scenario added (for adjusted prevalence estimation situation).
* Ability to specify number of optimisation iterations for `identifyBetaPars`
* `identifyBetaPars` is no longer an exported function; strictly internal
* Added support for normal and Poisson distributions in addition the beta distribution
* Added `getBetaFromCI`, a function that returns beta distribution functions from a given confidence interval (and in particular a random sampling function)
* Added `getNormFromCI`, a function that returns normal distribution functions from a given confidence interval (and in particular a random sampling function)
* Added `getPoisFromCI`, a function that returns Poisson distribution functions from a given confidence interval (and in particular a random sampling function)
* Changed the organisation of code within the source R code files
* Implemented changes requested from the submission of v0.1.0 to CRAN (references in DESCRIPTION file, TRUE/FALSE instead of T/F, restoring user parameters after doing plots)

# bootComb 0.1.0

* First version to be submitted to CRAN (but rejected).
