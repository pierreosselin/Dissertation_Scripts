source("InferenceCovariates.R")
source("Utils.R")
source("GraphModule.R")
library("mcmc")

n <- 10
nbIterations <- 1000000
step <- 10
rho <- c(3,1,1,1,1,1,1,1,1,1)
lK <- c(1,1,1,1,1,1,1,1,1,1)
a <- runif(n, min = -1, max = 1)
b <- runif(n, min = -1, max = 1)
sig <- 0
r <- 0.5
lambda <- rnorm(n*(n-1), mean = 0, sd = 1)
print(lambda)
vectParam <- c(a, b, sig, r, lambda)

chain <- run_metropolis_MCMC(vectParam, nbIterations, rho, lK, step)
burnIn = 10000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

par(mfrow = c(2,4))
hist(chain[-(1:burnIn),1],nclass=30,main="Sociability of 1", xlab="Mean value = red line" )
abline(v = mean(chain[-(1:burnIn),1]), col = "red")
hist(chain[-(1:burnIn),2],nclass=30, main="Sociability of 2", xlab="Mean value = red line")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")
hist(chain[-(1:burnIn),n + 1],nclass=30, main="Popularity of 1", xlab="Mean value = red line")
abline(v = mean(chain[-(1:burnIn),n+1]), col = "red")
hist(chain[-(1:burnIn),n + 2],nclass=30, main="Popularity of 2", xlab="Mean value = red line")
abline(v = mean(chain[-(1:burnIn),n+2]), col = "red")

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a")
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b")
plot(chain[-(1:burnIn),n+1], type = "l", xlab="True value = red line" , main = "Chain values of sd")
plot(chain[-(1:burnIn),n+2], type = "l", xlab="True value = red line" , main = "Chain values of sd")
# par(mfrow = c(1,1))
# visualizeSimpleGraph(rho, lK)
