## Test POsterior EM
source("Utils.R")
source("EM_Algorithm.R")
source("EM_VariableK.R")
n <- 20
K <- 5
p <- 4
a <- 2
b <- 2
lK <- rep(c(K), times = n)
rho <- generateNetwork(n,a,b,p,K)
rhoV <- as.vector(t(rho))
print(rho)
print(rhoV)
w <- matrix(rgamma(n*p, a, rate = b), n, p)
fit <- EMInference(rho,p,a,b, W = w, tol = 0.1)
fitV <- EMInferenceV(rhoV, lK, p, a, b, W = w, tol = 0.1)
lPost <- fit[["postList"]]
lPostV <- fitV[["postList"]]
plot((1:length(lPost))/100, lPost, type="l", xlab="Number of epochs", ylab="Log Posterior")
lines((1:length(lPostV))/100, lPostV, type="l", xlab="Number of epochs", ylab="Log Posterior", col = "red")

for (i in 1:20){
  fit <- EMInference(rho,p,a,b, tol = 0.1)
  lPost <- fit[["postList"]]
  lines((1:length(lPost))/100, lPost)
}
for (i in 1:20){
  fit <- EMInferenceV(rhoV, lK,p,a,b, tol = 0.1)
  lPost <- fit[["postList"]]
  lines((1:length(lPost))/100, lPost,  col = "red")
}
