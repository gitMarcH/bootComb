# bootComb 0.2.0.9000 (current development version)

* Added cranlogs badges to README.md
* Changed the way parameters are found for the poisson distribution. Rather than quantile matching, the CDF is matched (as this will be a continuous function rather than a step function in the way the quantile function is).
* Added support for negative binomial, gamma and exonential distributions in addition to the beta, normal and Poisson distributions. Finding the best-fit negative binomial distribution proceeds similar to the procedure for the Poisson distribution as this is another discrete distribution and hence the quantile function is a step function (albeit a two-dimensional step function in this case).
* Renamed the (internal only) function 'ssPoisPar' to 'ssPoisPars' for consistency with other distributions.
* Ranmed arguments 'pLow' and 'pUpp' from getBetaFromCI() to 'qLow' and 'qUpp' for consistency. Updated code from other functions where these parameters were specified.

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
