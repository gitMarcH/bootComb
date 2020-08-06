#' @title Compute the sum of squares between the theoretical and observed quantiles of a beta distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a beta distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit beta distribution for a given confidence interval.
#'
#' @param abPars The shape1 and shape2 parameters of the theoretical beta distribution.
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{qbeta}}
#'

ssBetaPars<-function(abPars,pLow,pUpp,alpha=0.05){
  res<-(qbeta(alpha/2,abPars[1],abPars[2])-pLow)^2+(qbeta(1-alpha/2,abPars[1],abPars[2])-pUpp)^2
  return(res)
}

#' @title Determine the best-fit beta distribution for a given confidence interval for a probability parameter.
#'
#' @description
#' Finds the best-fit beta distribution for a given confidence interval for a probability parameter and returns the shape1
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#'
#' @return A vector of length 2 giving the 2 parameter shape1 and shape1 for use with rebta/dbeta/ppbeta/qbeta
#'
#' @seealso
#' \code{\link{ssBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'
#' @examples
#' identifyBetaPars(pLow=0.1167,pUpp=0.1636,initPars=c(200,800))
#'
#' @export identifyBetaPars

identifyBetaPars<-function(pLow,pUpp,alpha=0.05,initPars=c(50,50)){
  if(pLow<0 | pUpp>1 | pLow>pUpp){stop("pLow and pUpp need to be both contained within [0,1] and pLow needs to be lower than pUpp.")}

  res<-suppressWarnings(optim(fn=ssBetaPars,pLow=pLow,pUpp=pUpp,par=initPars))
  if(res$convergence!=0){stop("optim() as called by identifyBetaPars() failed to converge.")}

  return(res$par)
}

#' @title Simulation scenario 1.
#'
#' @description
#' This is a simulation to compute the coverage of the confidence interval returned by bootComb() in the case of the product of 2 probability parameter estimates.
#'
#' @param B The number of simulations to run. Defaults to 1e3.
#' @param p1 The true value of the first probability parameter.
#' @param p2 The true value of the second probability parameter.
#' @param nExp1 The size of each simulated experiment to estimate p1.
#' @param nExp2 TThe size of each simulated experiment to estimate p2.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the proportion of simulations for which the confidence interval contained the true parameter value.
#'
#' @examples
#' \dontrun{
#' simScen1(p1=0.35,p2=0.2,nExp1=100,nExp2=1000,B=100)
#'   # B value only for convenience here
#'   # Increase B to 1e3 or 1e4 (be aware this may run for some time).
#'  }
#'
#' @export simScen1

simScen1<-function(B=1e3,p1,p2,nExp1,nExp2,alpha=0.05){
  trueP<-p1*p2
  res<-rep(0,B)

  for(j in 1:B){
    exp1<-rbinom(1,size=nExp1,prob=p1)
    exp2<-rbinom(1,size=nExp2,prob=p2)

    p1Est<-binom.test(x=exp1,n=nExp1)
    p2Est<-binom.test(x=exp2,n=nExp2)

    betaPars1<-identifyBetaPars(pLow=p1Est$conf.int[1],pUpp=p1Est$conf.int[2],initPars=c(round(nExp1/2),round(nExp1/2)))
    betaPars2<-identifyBetaPars(pLow=p2Est$conf.int[1],pUpp=p2Est$conf.int[2],initPars=c(round(nExp2/2),round(nExp2/2)))

    dist1<-function(n){rbeta(n,betaPars1[1],betaPars1[2])}
    dist2<-function(n){rbeta(n,betaPars2[1],betaPars2[2])}

    distListEx<-list(dist1,dist2)
    combFunEx<-function(pars){pars[[1]]*pars[[2]]}

    bsOut<-bootComb(distList=distListEx,combFun=combFunEx)

    if(trueP>=bsOut$conf.int[1] & trueP<=bsOut$conf.int[2]){res[j]<-1}
  }

  coverage<-sum(res)/length(res) # ideally should be 95%

  return(coverage)
}
