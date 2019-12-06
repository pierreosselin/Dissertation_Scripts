## Compute the mean difference in ranking
rankingGap <- function(v1,v2){
  n <- length(v1)
  l1 <- 1:n
  l2 <- match(v1,v2)
  l3 <- l1 - l2
  l3 <- abs(l3)
  res <- mean(l3)
  return(res)
}

## Generate list Plackett-Luce Model
listPlackett <- function(n, K, lambda){
  rho <- c()
  for (i in 1:K){
    prob <- lambda / sum(lambda)
    res <- rmultinom(1,1, prob = prob)
    ind <- match(c(1), res)
    lambda[ind] <- 0
    rho <- c(rho, ind)
  }
  return(rho)
}

listPlackettL <- function(n, K, lambda){
  rho <- c()
  for (i in 1:K){
    prob <- lambda / sum(lambda)
    res <- rmultinom(1,1, prob = prob)
    ind <- match(c(1), res)
    lambda[ind] <- 0
    rho <- c(rho, ind)
  }
  return(rho)
}

## Generate Networks with community representation

generateNetwork <- function(n,a,b,p,K){
  W <- matrix(rgamma(n*p, a, rate = b), n, p)
  lambda <- W %*% t(W)
  for (i in 1:n){
    lambda[i,i] <- 0
  }
  rho <- apply(lambda, 1, (function(x) listPlackett(n, K, x))) 
  return(t(rho))
}

generateNetworkL <- function(n,a,b, listK){
  lambda <- matrix(rgamma(n*n, a, rate = b), n, n)
  for (i in 1:n){
    lambda[i,i] <- 0
  }
  rho <- c()
  for (i in 1:n){
    rho <-c(rho, listPlackett(n, listK[i], lambda[i,]))
  }
  return(rho)
}

plotResult <- function(res) {
  lPost <- res[["postList"]]
  plot((1:length(lPost))/100, lPost, type="l", xlab="Number of epochs", ylab="Log Posterior")
}