library("MASS", lib.loc="/usr/lib/R/library")
library("amen", lib.loc="/usr/lib/R/library")
sim <- 10
n <- 100
maxRank <- 100
Sigma_ab <- matrix(c(1,0.5,0.5,1),2,2)
Sigma_eps <- matrix(c(1,0.9,0.9,1),2,2)
sample_ab <- mvrnorm(n = n, rep(0, 2), Sigma_ab)
sample_eps <- mvrnorm(n = n*(n-1)/(2), rep(0, 2), Sigma_eps)
beta <- 1
Y <- matrix(rep(0,n*n), n, n)
for (i in 1:n){
  for (j in 1:n){
    if(i == j){
      Y[i,i] <- NA
    } else if (i < j){
      Y[i,j] <- beta + sample_ab[i] + sample_ab[j + n] + sample_eps[(i-1)*(n - (i/2)) + j - i]
    } else {
      Y[i,j] <- beta + sample_ab[i] + sample_ab[j + n] + sample_eps[(j-1)*(n - (j/2)) + i - j + n*(n-1)/(2)]
    }
  }
}
Yrank <- matrix(rep(0,n*n), n, n)
for (i in 1:n){
  Yrank[i,i] <- NA
  lRank <- order(Y[i, ], na.last = NA, decreasing = TRUE)
  lRank <- lRank[Y[i, lRank[]] > 0]
  indice <- min(maxRank, length(lRank))
  for (j in 1:indice){
    Yrank[i,lRank[j]] <- indice - j + 1
  }
}

fit<-ame(Y,R=2,model = "frn",odmax=4)

summary(fit)
