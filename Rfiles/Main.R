## Main Script
source("Utils.R")
source("EM_Algorithm.R")
data("sampson")
n <- 20
K <- 5
p <- 4
a <- 2
b <- 2
rho <- generateNetwork(n,a,b,p,n-1)
print(rho)
Wresult <- EMInference(rho[,1:K],p,a,b, tol = 0.000001)
print(Wresult)

lambdainf <- Wresult %*% t(Wresult)
for (i in 1:n){
  lambdainf[i,i] <- 0
}
tes <- listPlackett(n,n-1,lambdainf[1,])
print(rankingGap(tes, rho[1,]))
