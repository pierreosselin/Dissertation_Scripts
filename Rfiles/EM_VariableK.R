# R Script performing EM algorithm for Community Representation, in this algorithm the K is not fixed.
###


## Posterior Likelihood of the (W)
PosteriorLikeV <- function(W, a, b, rho, lK){
  n <- length(lK)
  shape <- dim(W)
  p <- shape[2]
  remove(shape)
  lambda <- W %*% t(W)
  for (i in 1:n){
    lambda[i,i] <- 0
  }
  
  Delta <- rep(rowSums(lambda), times = lK)
  for (i in 1:n){
    if (i == 1) {
      if (lK[i] > 1){
        for (j in 2:lK[i]){
          Delta[j] <- Delta[j-1] - lambda[i, rho[j-1]]
        }
      }
    } else {
      if (lK[i] > 1){
        for (j in 2:lK[i]){
          Delta[sum(lK[1:(i-1)]) + j] <- Delta[sum(lK[1:(i-1)]) + j - 1] - lambda[i, rho[sum(lK[1:(i-1)]) + j - 1]]
        }
      }
    }
  }
  
  interm <- c()
  for (j in 1:(lK[1])){
    if (lK[i] > 0) {
      interm <- c(interm, lambda[1, rho[j]])
    }
  }
  for (i in 2:n){
    if (lK[i] > 0) {
      for(j in 1:(lK[i])){
        interm <- c(interm, lambda[i,rho[(sum(lK[1:(i-1)]) + j)]])
      }
    }
  }
  
  L <- sum( (a-1) * (log(W)) - b * W ) - sum(log(Delta)) + sum(log(interm))
  return(L)
}


## EM algorithm with variable rank, lK is the list of number of ranks for each individual
EMInferenceV <- function(rho,lK, p, a, b, W = FALSE, tol = 0.1) {
  lineR <- function(r) {
    if (r == 1) {
      indexLine <- 1:lK[1]
    } else {
      indexLine <- (sum(lK[1:(r-1)]) + 1):(sum(lK[1:r]))
    }
    return(indexLine)
  }
  
  # Update function for w
  update_rs <- function(r, s, rho, lK, n, p, W, a, b, lambda, invDelta) {
    #Compute C
    
    
    l <- 1:n
    l <- l[-r]
    indexLineR <- lineR(r)
    if (lK[r] > 0){
      C <- a - 1 + W[r,s] * sum(W[rho[indexLineR], s] /  lambda[r,rho[indexLineR]]) +  W[r,s] * ((W[l,s] / lambda[l,r] ) %*% (sapply(l, (function (x)   if (lK[x] > 0) {    (if (x == 1) {as.numeric(r %in% rho[1:lK[1]])} else {as.numeric(r %in% rho[(sum(lK[1:(x-1)]) + 1):(sum(lK[1:x]))])})    } else {  0  }        ) ) ))
      A1 <- b + (sum(invDelta[indexLineR])) * (sum(W[l,s]))
      A2 <- (W[l,s] %*% sapply(l, (function(x)  if (lK[x] > 0)  { (sum(invDelta[lineR(x)]) )  } else {0}       ) ))
      A3 <- if (lK[r] > 1)  {(-1)*(W[rho[head(indexLineR, -1)] ,s] %*% sapply(1:(lK[r] - 1), (function(i) sum(invDelta[if (r == 1) {(1 + i):lK[1]} else {(sum(lK[1:(r-1)]) + 1 + i):(sum(lK[1:r]))}]))))} else {0}
      A4 <- (-1)*(W[l,s] %*% sapply(l, (function(i)  if (lK[i] > 0) {   if (r %in% head(rho[lineR(i)], -1)) {(sum(invDelta[if (i == 1) {(match(r, rho[lineR(i)]) + 1):lK[1]} else {(sum(lK[1:(i-1)]) + 1 + match(r, rho[lineR(i)])):(sum(lK[1:i]))}]))} else {0}     } else {0}           )))
      A <- A1 + A2 + A3 + A4
    } else {
      C <- a - 1 + W[r,s] * ((W[l,s] / lambda[l,r] ) %*% (sapply(l, (function (x)     if (lK[x] > 0)  {   (if (x == 1) {as.numeric(r %in% rho[1:lK[1]])} else {as.numeric(r %in% rho[(sum(lK[1:(x-1)]) + 1):(sum(lK[1:x]))])})     }   else {0}       ) ) ))
      A <- b + (W[l,s] %*% sapply(l, (function(x)   if (lK[x] > 0)  { (sum(invDelta[lineR(x)]) )  } else {0} ) )) - (W[l,s] %*% sapply(l, (function(i)    if (lK[i] > 0)  {     if (r %in% head(rho[lineR(i)], -1)) {(sum(invDelta[if (i == 1) {(match(r, rho[lineR(i)]) + 1):lK[1]} else {(sum(lK[1:(i-1)]) + 1 + match(r, rho[lineR(i)])):(sum(lK[1:i]))}]))} else {0}   } else {0}      )))
    }
    newVal <- C / A
    change_rs <- abs(W[r,s]  - newVal)
    W[r,s] <- newVal
    return(list("W" = W, "amplitude" = change_rs))
  }
  
  # Compute basic quantities of interest
  n <- length(lK)
  totK <- sum(lK)
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
        Delta <- rep(rowSums(lambda), times = lK)
        for (i in 1:n){
          if (i == 1) {
            if (lK[i] > 1) {
              for (j in 2:lK[i]){
                Delta[j] <- Delta[j-1] - lambda[i, rho[j-1]]
              }
            }
          } else {
            if (lK[i] > 1){
              for (j in 2:lK[i]){
                Delta[sum(lK[1:(i-1)]) + j] <- Delta[sum(lK[1:(i-1)]) + j - 1] - lambda[i, rho[sum(lK[1:(i-1)]) + j - 1]]
              }
            }
          }
        }
        invDelta = 1 / Delta
        res <- update_rs(r, s, rho, lK, n, p, W, a, b, lambda, invDelta)
        maxCurrentEpoch <- max(maxCurrentEpoch, res[["amplitude"]])
        W <- res[["W"]]
        lPost <- c(lPost, PosteriorLikeV(W, a, b, rho, lK))
      }
    }
  }
  return(list("W" = W, "postList" = lPost))
}