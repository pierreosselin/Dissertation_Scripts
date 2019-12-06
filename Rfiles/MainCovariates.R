source("InferenceCovariates.R")
source("Rfiles/Utils.R")
library("mcmc")

n <- 6
K <- 3
p <- 2
a <- 2
b <- 2
nbIterations <- 50000
rho <- generateNetwork(n,a,b,p,K)
a <- runif(n, min = -1, max = 1)
b <- runif(n, min = -1, max = 1)
sig <- 0
r <- 0.5
lambda <- rnorm(n*(n-1), mean = 0, sd = 1)
vectParam <- c(a, b, sig, r, lambda)
print(vectParam)

lK <- rep(K, n)
rho <- as.vector(t(rho))
chain <- run_metropolis_MCMC(vectParam, nbIterations, rho, lK)
targetdistrib(vectParam, rho = rho, lK = lK)
burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=30,main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]), col = "red")
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), col = "red")
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a")
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b")
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd")
