source("InferenceCovariates.R")
source("Utils.R")
source("GraphModule.R")
library("mcmc")
library("lda")

data("sampson")

df <- sampson[["SAMPLK1"]]
n <- dim(df)[1]

mat <- matrix(df, n, n)
Kmax <- max(mat)
listK <- c()
rhoK <- c()

for (i in 1:n) {
  interm <- c()
  for (j in 1:Kmax) {
    interm <- c(interm, which(mat[i, ] == j))
  }
  rhoK <- c(rhoK, interm)
  listK <- c(listK, length(interm))
}


nbIterations <- 50000
a <- runif(n, min = -1, max = 1)
b <- runif(n, min = -1, max = 1)
sig <- 0
r <- 0.5
lambda <- rnorm(n*(n-1), mean = 0, sd = 1)
vectParam <- c(a, b, sig, r, lambda)


chain <- run_metropolis_MCMC(vectParam, nbIterations, rho = rhoK, lK = listK)
burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
matRes <- matrix(rep(0,2*n), n, 2)
for (i in 1:n){
  matRes[i,1] <- mean(chain[-(1:burnIn),i])
  matRes[i,2] <- mean(chain[-(1:burnIn),n+i])
}

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),19],nclass=30,main="Popularity of 1", xlab="Mean value = red line" )
abline(v = mean(chain[-(1:burnIn),19]), col = "red")
hist(chain[-(1:burnIn),31],nclass=30, main="Popularity of 13", xlab="Mean value = red line")
abline(v = mean(chain[-(1:burnIn),31]), col = "red")
hist(chain[-(1:burnIn),33],nclass=30, main="Popularity of 15", xlab="Mean value = red line")
abline(v = mean(chain[-(1:burnIn),33]), col = "red")

plot(chain[-(1:burnIn),19], type = "l", xlab="True value = red line" , main = "Chain values of a")
plot(chain[-(1:burnIn),31], type = "l", xlab="True value = red line" , main = "Chain values of b")
plot(chain[-(1:burnIn),33], type = "l", xlab="True value = red line" , main = "Chain values of sd")
par(mfrow = c(1,1))
visualizeSimpleGraph(rhoK, listK)
