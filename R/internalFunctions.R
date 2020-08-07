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
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
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

identifyBetaPars<-function(pLow,pUpp,alpha=0.05,initPars=c(50,50),maxiter=1e3){
  if(pLow<0 | pUpp>1 | pLow>pUpp){stop("pLow and pUpp need to be both contained within [0,1] and pLow needs to be lower than pUpp.")}

  res<-suppressWarnings(optim(fn=ssBetaPars,pLow=pLow,pUpp=pUpp,par=initPars,control=list(maxit=maxiter)))
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
#' @param nExp1 The size of each simulated experiment to estimate \code{p1}.
#' @param nExp2 TThe size of each simulated experiment to estimate \code{p2}.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A list with 2 elements:
#'    \code{estimate} A single number, the proportion of simulations for which the confidence interval contained the true parameter value.
#'    \code{conf.int} A 95% confidence interval for the coverage estimate.
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

  tmp<-binom.test(x=sum(res),n=length(res))
  coverage<-list(estimate=tmp$estimate,conf.int=tmp$conf.int) # ideally should be 95%

  return(coverage)
}

#' @title Simulation scenario 2.
#'
#' @description
#' This is a simulation to compute the coverage of the confidence interval returned by bootComb() in the case of adjusting a prevalence estimate for estimates of sensitivity and specificity.
#'
#' @param B The number of simulations to run. Defaults to 1e3.
#' @param p The true value of the prevalence parameter.
#' @param sens The true value of the assay sensitivity parameter.
#' @param spec The true value of the assay specificity parameter
#' @param nExp The size of each simulated experiment to estimate \code{p}.
#' @param nExpSens The size of each simulated experiment to estimate \code{sens}.
#' @param nExpSpec The size of each simulated experiment to estimate \code{spec}.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A list with 2 elements:
#'    \code{estimate} A single number, the proportion of simulations for which the confidence interval contained the true prevalence parameter value.
#'    \code{conf.int} A 95% confidence interval for the coverage estimate.
#'
#' @examples
#' \dontrun{
#' simScen2(p=0.15,sens=0.90,spec=0.95,nExp=250,nExpSens=1000,nExpSpec=500,B=100)
#'   # B value only for convenience here
#'   # Increase B to 1e3 or 1e4 (be aware this may run for some time).
#'  }
#'
#' @export simScen2

simScen2<-function(B=1e3,p,sens,spec,nExp,nExpSens,nExpSpec,alpha=0.05){
  trueObsPrev<-p*sens+(1-p)*(1-spec)
  res<-rep(0,B)

  for(j in 1:B){
    expP<-rbinom(1,size=nExp,prob=trueObsPrev)
    expSens<-rbinom(1,size=nExpSens,prob=sens)
    expSpec<-rbinom(1,size=nExpSpec,prob=spec)

    pEst<-binom.test(x=expP,n=nExp)
    sensEst<-binom.test(x=expSens,n=nExpSens)
    specEst<-binom.test(x=expSpec,n=nExpSpec)

    bsOut<-adjPrevSensSpecCI(prevCI=pEst$conf.int,sensCI=sensEst$conf.int,specCI=specEst$conf.int,alpha=alpha)

    if(p>=bsOut$conf.int[1] & p<=bsOut$conf.int[2]){res[j]<-1}
  }

  tmp<-binom.test(x=sum(res),n=length(res))
  coverage<-list(estimate=tmp$estimate,conf.int=tmp$conf.int) # ideally should be 95%

  return(coverage)
}
