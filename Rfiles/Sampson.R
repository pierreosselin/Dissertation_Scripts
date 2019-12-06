## Main Script
source("Utils.R")
source("EM_Algorithm.R")
source("EM_VariableK.R")
source("GraphModule.R")
library("lda")
data(sampson)
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
print(rhoK)
print(listK)
p <- 3
a <- 2
b <- 2
Wresult1 <- EMInferenceV(rhoK, listK,p,a,b, tol = 0.01)
# plotResult(Wresult)
# visualizeGraph(Wresult, rhoK, listK)


df <- sampson[["SAMPLK2"]]
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
print(rhoK)
print(listK)
p <- 3
a <- 2
b <- 2
Wresult2 <- EMInferenceV(rhoK, listK,p,a,b, tol = 0.01)
plotResult(Wresult)
visualizeGraph(Wresult, rhoK, listK)



df <- sampson[["SAMPLK3"]]
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
print(rhoK)
print(listK)
p <- 3
a <- 2
b <- 2
Wresult3 <- EMInferenceV(rhoK, listK,p,a,b, tol = 0.01)
plotResult(Wresult)
visualizeGraph(Wresult, rhoK, listK)

lPost1 <- Wresult1[["postList"]]
lPost2 <- Wresult2[["postList"]]
lPost3 <- Wresult3[["postList"]]

plot((1:length(lPost1))/100, lPost1, type="l", xlab="Number of epochs", ylab="Log Posterior", col = "blue")
lines((1:length(lPost2))/100, lPost2, type="l", xlab="Number of epochs", ylab="Log Posterior", col = "red")
lines((1:length(lPost3))/100, lPost3, type="l", xlab="Number of epochs", ylab="Log Posterior", col = "green")
legend(1, 95, legend=c("Line 1", "Line 2", "Line 3"), col=c("blue", "red", "green"), lty=1:2, cex=0.8)
