# R Script performing EM algorithm for Community Representation
###


## Posterior Likelihood of the (W)
PosteriorLike <- function(W, a, b, rho){
  shape <- dim(rho)
  n <- shape[1]
  K <- shape[2]
  shape <- dim(W)
  p <- shape[2]
  remove(shape)
  
  lambda <- W %*% t(W)
  for (i in 1:n){
    lambda[i,i] <- 0
  }
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
  L <- sum( (a-1) * (log(W)) - b * W ) - sum(log(Delta)) + sum(log(interm))
  return(L)
}


## EM algorithm
EMInference <- function(rho, p, a, b, W = FALSE, tol = 0.1) {
  
  # Update function for w
  update_rs <- function(r, s, rho, n, p, K, W, a, b, lambda, invDelta) {
    #Compute C
    l <- 1:n
    l <- l[-r]
    C <- a - 1 + W[r,s] * sum(W[rho[r,], s] /  lambda[r,rho[r,]]) +  W[r,s] * (W[l,s] %*% (apply(rho[l,], 1, (function(x) as.integer(r %in% x)) ) / lambda[l,r]  ))
    A <- b + (sum(invDelta[r,])) * (sum(W[l,s])) + (W[l,s] %*% rowSums(invDelta[l,])) - (W[rho[r, 1:(K-1)] ,s] %*% sapply(1:(K-1), (function(i) sum(invDelta[r, (i+1):K])))) - (W[l,s] %*% sapply(l, (function(i) if (r %in% rho[i, 1:(K-1)]) (sum(invDelta[i, (match(r, rho[i, ]) + 1):K])) else 0 )))
    newVal <- C / A
    change_rs <- abs(W[r,s]  - newVal)
    W[r,s] <- newVal
    return(list("W" = W, "amplitude" = change_rs))
  }
  
  # Compute basic quantities of interest
  shape <- dim(rho)
  n <- shape[1]
  K <- shape[2]
  remove(shape)

  lPost <- c()
  maxCurrentEpoch <- tol + tol/1000
  security <- 0
  maxSecurity <- 5000
  
  if (is.numeric(W) == FALSE){
    W <- matrix(rgamma(n*p, a, rate = b), n, p)
  }
  
  while (maxCurrentEpoch > tol & security < maxSecurity)
  {
    maxCurrentEpoch <- 0
    security <- security + 1
    # Perform the whole W update
    for (r in 1:n){
      for (s in 1:p){
        lambda <- W %*% t(W)
        for (i in 1:n){
          lambda[i,i] <- 0
        }
        Delta <- matrix(rowSums(lambda), n, K)
        for (i in 1:n){
          for (j in 2:K){
            Delta[i,j] <- Delta[i,j-1] - lambda[i, rho[i,j-1]]
          }
        }
        invDelta <- 1 / Delta
        res <- update_rs(r, s, rho, n, p, K, W, a, b, lambda, invDelta)
        maxCurrentEpoch <- max(maxCurrentEpoch, res[["amplitude"]])
        W <- res[["W"]]
        lPost <- c(lPost, PosteriorLike(W, a, b, rho))
      }
    }
    print(maxCurrentEpoch)
  }
  return(list("W" = W, "postList" = lPost))
}