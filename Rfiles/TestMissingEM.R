## Script for testing missing data in EM algo

source("Utils.R")
source("EM_Algorithm.R")
source("EM_VariableK.R")
source("GraphModule.R")

n <- 6
p <- 2
a <- 2
b <- 2
lK <- c(1,1,1,1,1,0)
rho <- c(2,1,1,1,1,1,1)
print(rho)
w <- matrix(rgamma(n*p, a, rate = b), n, p)
fitV <- EMInferenceV(rho, lK, p, a, b, W = w, tol = 0.1)
lPostV <- fitV[["postList"]]
wRes <- fitV[["W"]]
print(wRes)
plot((1:length(lPostV))/12, lPostV, type="l", xlab="Number of epochs", ylab="Log Posterior")
#lines((1:length(lPostV))/100, lPostV, type="l", xlab="Number of epochs", ylab="Log Posterior", col = "red")
lambda <- wRes %*% t(wRes)
for (i in 1:n){
  lambda[i,i] <- 0
}
print(lambda)

#visualizeSimpleGraph(rho, lK)
