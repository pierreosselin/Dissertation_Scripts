library("MASS", lib.loc="/usr/lib/R/library")
library("mvtnorm", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("invgamma")
library("progress")
## First, there will be no covariates beta = 0

# s = prior in beta, a, b. r = prior in the inverse wishart
loglikelihood <- function(lambda, rho) {
  shape <- dim(rho)
  n <- shape[1]
  K <- shape[2]
  remove(shape)
  
  Delta <- matrix(rowSums(lambda), n, K)
  for (i in 1:n){
    for (j in 2:K){
      Delta[i,j] <- Delta[i,j-1] - lambda[i, rho[i,j-1]]
    }
  }
  interm <- matrix(0,n,K)
  for (i in 1:n){
    for(j in 1:K){
      interm[i,j] <- lambda[i,rho[i,j]]
    }
    }
  L <- sum(sum(log(Delta)) + sum(log(interm)))
}
dataFunction <- function(){
  n <- 10
  priorParam <- 10**(5)
  priorSigmaAlpha <- 1
  priorSigmaBeta <-0.01
  lengthCov <- 0
  return(list("n" = n, "priorParam" = priorParam, "lengthCov" = lengthCov, "priorSigmaAlpha" = priorSigmaAlpha, "priorSigmaBeta" = priorSigmaBeta))
}

#vecParam = [a, b, sigma, rho, lambda, beta] Here lambda are log(lambda) and without the diagonal and sigma and rho are the log of the real parameters
targetdistrib <- function(vecParam, rho = rho, lK = lK) {
  dat <- dataFunction()
  n <- dat[["n"]]
  pr <- dat[["priorParam"]]
  prSigmaAlpha <- dat[["priorSigmaAlpha"]]
  prSigmaBeta <- dat[["priorSigmaBeta"]]
  lengthCov <- dat[["lengthCov"]]
  a <- vecParam[1:n]
  b <- vecParam[(n+1):(2*n)]
  sig <- vecParam[(2*n + 1)]
  realr <- vecParam[(2*n + 2)]
  realSigma <- exp(sig)
  lambdavec <- matrix(vecParam[(2*n + 3):(2*n + 2 + n*(n-1))], n, (n-1), byrow = TRUE)
  interm <- matrix(0,n,n)
  for (i in 1:n){
    interm[i,-i] <- lambdavec[i,]
  }
  lambda <- interm
  condBeta <- 0
  if (lengthCov != 0) {
    beta <- vecParam[-1:-(2*n + 3 + n*(n-1))]
    condBeta <- -sum((beta**2)/(2*(pr**2)))
  }
  
  #Computation prio
  prior <- -sum((a**2)/(2*(pr**2))) - sum((b**2)/(2*(pr**2))) + condBeta - (prSigmaAlpha - 1)*sig - (prSigmaBeta)*realSigma
  
  #Computation snd term
  summation <- 0
  for (i in (1:(n-1))) {
    for (j in (i+1):n) {
      summation <- summation + dmvnorm(c(lambda[i,j], lambda[j,i]), mean = c(a[i] + b[j], a[j] + b[i]),sigma= matrix(c(realSigma**2, (realSigma**2)*realr, (realSigma**2)*realr, realSigma**2), 2,2), log=TRUE)
    }
  }
  
  
  #Computation 3nd term
  lambdaExp <- exp(lambda)
  Delta <- rep(rowSums(lambdaExp), times = lK)
  for (i in 1:n){
    if (i == 1) {
      if(lK[i] > 1){
        for (j in 2:lK[i]){
          Delta[j] <- Delta[j-1] - lambdaExp[i, rho[j-1]]
        }
      }
    } else {
      if (lK[i] > 1){
        for (j in 2:lK[i]){
          Delta[sum(lK[1:(i-1)]) + j] <- Delta[sum(lK[1:(i-1)]) + j - 1] - lambdaExp[i, rho[sum(lK[1:(i-1)]) + j - 1]]
        }
      }
    }
  }
  
  interm <- c()
  for (j in 1:(lK[1])){
    interm <- c(interm, lambdaExp[1, rho[j]])
  }
  for (i in 2:n){
    for(j in 1:(lK[i])){
      interm <- c(interm, lambdaExp[i,rho[(sum(lK[1:(i-1)]) + j)]])
    }
  }
  
  
  
  
  
  # Delta <- matrix(rowSums(lambdaExp), n, K)
  # for (i in 1:n){
  #   for (j in 2:K){
  #     Delta[i,j] <- Delta[i,j-1] - lambdaExp[i, rho[i,j-1]]
  #   }
  # }
  # 
  # interm <- matrix(0,n,K)
  # for (i in 1:n){
  #   for(j in 1:K){
  #     interm[i,j] <- lambda[i,rho[i,j]]
  #   }
  # }
  
  thirdTerm <- -sum(log(Delta)) + sum(log(interm))
  return(prior + summation + thirdTerm)
}


proposalfunction <- function(vecParam){
  dat <- dataFunction()
  n <- dat[["n"]]
  lengthCov <- dat[["lengthCov"]]
  if (runif(1) < 0.1){
    l1 = rnorm(2*n + 1, mean = rep(0, 2*n + 1), sd = rep(5, 2*n + 1))
    l2 = runif(1)
    l3 = rnorm(n*(n-1), mean = rep(0, n*(n-1)), sd = rep(5,n*(n-1)))
    return(c(l1,l2,l3))
  }
  r = vecParam[(2*n + 2)]
  lambdavec <- matrix(vecParam[(2*n + 3):(2*n + 2 + n*(n-1))], n, (n-1), byrow = TRUE)
  l1 = rnorm(2*n + 1, mean = vecParam[1:(2*n + 1)], sd = rep(0.125, 2*n + 1))
  l2 = runif(1, min = max(r - 0.075, 0), max = min(1, r + 0.075))
  l3 = rnorm(n*(n-1), mean = vecParam[(2*n + 3):(2*n + 2 + n*(n-1))], sd = rep(0.125,n*(n-1)))
  return(c(l1,l2,l3))
}

run_metropolis_MCMC <- function(startvalue, iterations, rho, lK, step = 1){
  pb <- progress_bar$new(total = iterations)
  d <- length(startvalue)
  chain = array(dim = c((iterations/10)+1,d))
  chain[1,] = startvalue
  currentPoint = startvalue
  for (i in 1:iterations){
    #proposal = proposalfunction(chain[i,])
    #probab = exp(targetdistrib(proposal, rho = rho, lK = lK) - targetdistrib(chain[i,], rho = rho, lK = lK))
    proposal = proposalfunction(currentPoint)
    probab = exp(targetdistrib(proposal, rho = rho, lK = lK) - targetdistrib(currentPoint, rho = rho, lK = lK))
    if (runif(1) < probab){
      #chain[i+1,] = proposal
      currentPoint <- proposal
    }
    #else{
      #chain[i+1,] = chain[i,]
    #}
    if (i %% step == 0){
      chain[(i %/% step) + 1,] = proposal
    }
    pb$tick()
  }
  return(chain)
}