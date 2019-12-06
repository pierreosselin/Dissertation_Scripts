## Test POsterior EM
source("Utils.R")
source("EM_Algorithm.R")
source("EM_VariableK.R")
n <- 25
K <- 5
p <- 4
a <- 2
b <- 1
nSimu <- 10
nSampleRank <- 20
estimateError <- 0
for (l in 1:nSimu){
  lK <- rep(c(K), times = n)
  rho <- generateNetwork(n,a,b,p,K)
  rhoV <- as.vector(t(rho))
  w <- matrix(rgamma(n*p, a, rate = b), n, p)
  lambdaReal <- w %*% t(w)
  fitV <- EMInferenceV(rhoV, lK, p, a, b, W = w, tol = 0.1)
  wRes <- fitV[["W"]]
  lambda <- wRes %*% t(wRes)
  for (i in 1:n){
    lambda[i,i] <- 0
    lambdaReal[i,i] <- 0
  }
  error <- 0
  for (j in 1:nSampleRank){
    for (i in 1:n){
      inferedRank <- listPlackett(n,n-1,lambda[i,])
      realRank <- listPlackett(n,n-1,lambdaReal[i,])
      interm <- rankingGap(realRank, inferedRank)
      error <- error + interm
    }
  }
  estimateError <- estimateError + error/(n*nSampleRank)
  print("Simulation Done")
}
estimateError <- estimateError/nSimu
