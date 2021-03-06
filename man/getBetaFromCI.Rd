% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitDistributionsToCIs.R
\name{getBetaFromCI}
\alias{getBetaFromCI}
\title{Find the best-fit beta distribution for a given confidence interval for a probability parameter.}
\usage{
getBetaFromCI(qLow, qUpp, alpha = 0.05, initPars = c(50, 50), maxiter = 1000)
}
\arguments{
\item{qLow}{The observed lower quantile.}

\item{qUpp}{The observed upper quantile.}

\item{alpha}{The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.}

\item{initPars}{A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).}

\item{maxiter}{Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.}
}
\value{
A list with 5 elements:
\item{r}{The sampling function.}
\item{d}{The density function.}
\item{p}{The distribution function.}
\item{q}{The quantile function.}
\item{pars}{A vector of length 2 giving the two shape parameters for the best-fit beta distribution (\code{shape1} and \code{shape2} as in \code{\link{rbeta}}, \code{\link{dbeta}}, \code{\link{pbeta}}, \code{\link{qbeta}}).}
}
\description{
Finds the best-fit beta distribution for a given confidence interval for a probability parameter; returns the corresponding density, distribution, quantile and sampling functions.
}
\examples{
b<-getBetaFromCI(qLow=0.1167,qUpp=0.1636,initPars=c(200,800))
print(b$pars) # the fitted parameter values
b$r(10) # 10 random values from the fitted beta distribution
b$d(0.15) # the probability density at x=0.15 for the fitted beta distribution
b$p(0.15) # the cumulative density at x=0.15 for the fitted beta distribution
b$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
x<-seq(0,1,length=1e3)
y<-b$d(x)
plot(x,y,type="l",xlab="",ylab="density") # density plot for the fitted beta distribution

}
\seealso{
\code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
}
