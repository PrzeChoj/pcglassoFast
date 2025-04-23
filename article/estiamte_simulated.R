library(pcglassoFast)
library(glasso)
library(ggplot2)
library(Matrix)
#graphics.off()
generate.pcglasso <- TRUE #either pcpglasso or glasso
# Number of observations you want to simulate
n <- 1000000

# simulate data from
if(generate.pcglasso==FALSE){
  data(Q_simulated_glasso)
  Q <- Q_simulated_glasso
  lambda.pc <- 3 #0.9995
  lambda.glasso <- 0.05
}else{
  data(Q_simulated_pcglasso)
  Q <- Q_simulated_pcglasso
  lambda.pc <- 0.11
  lambda.glasso <- 0.025

}
L <- Cholesky(Matrix(forceSymmetric(Q), sparse = TRUE), LDL = FALSE, perm = TRUE)


dim <- ncol(Q)


# Simulate standard normal data
z <- matrix(rnorm(n * dim), nrow = dim, ncol = n)

# Solve Q^(1/2) x = z for x, which gives samples from N(0, Q^-1)
simulated_data <- Matrix::solve(L,  Matrix::solve(L, z, system = "P"), system = "Lt")
simulated_data <- Matrix::solve(L, simulated_data, system = "Pt")
# Transpose to have observations in rows
simulated_data_matrix <- as.matrix(t(simulated_data))
S <- cov(simulated_data_matrix)
lambdas <- c(0.12,0.08)
pcPath <- pcglassoPath(S,alpha=0, max.edge.fraction = 0.4)
pcPath_rev <- pcglassoPath(S,alpha=0,lambdas = rev(pcPath$lambdas), max.edge.fraction = 1)
plot(pcPath$lambdas,pcPath$loss,col='red',type='l')
lines(pcPath_rev$lambdas,pcPath_rev$loss,col='blue')
loss_pcPath <- loss.evaluation(pcPath, Sigma=S,n=n, gamma=0.)
Q_prun <- solve(S)
Q_prun[abs(Q_prun)< 1e-1]=0
image(Matrix(Q,sparse=T))

L.P <- chol(Q_prun)
p <- dim(Q_prun)[1]
gamma <- 0
loglik <- (sum(log(diag(L.P))) - 0.5 * sum(diag(Q_prun%*%S)) - 0.5 * p * log(2 * pi))

loglik_true <-  (sum(log(diag(chol(Q)))) - 0.5 * sum(diag(Q%*%S)) - 0.5 * p * log(2 * pi))
n_edges <- (p * p - sum(Q==0) - p)/2 + p
bic_true <- -2 * loglik_true +n_edges * (log(n) + 4 * gamma * log(p))
