#need package space installed through
#remotes::install_version("space", version = "0.1-1.1",
#                         repos = "https://cran.r-project.org")
library(space)
library(pcglassoFast)
source("article/simulation_functions.R")
library(glasso)
library(ggplot2)
library(Matrix)
set.seed(2)
graphics.off()
generate.pcglasso <- FALSE #either pcpglasso or glasso
split.train <- 0.7
# Number of observations you want to simulate
ns=c(200,300)
sim = 2
nlambda <- 100
mc_cores = 6
alpha.grid=0
lambda.min.ratio <- 0.01
# simulate data from
if(generate.pcglasso==FALSE){
  data(Q_simulated_glasso)
  Q <- Q_simulated_glasso
}else{
  data(Q_simulated_pcglasso)
  Q <- Q_simulated_pcglasso

}
Q <- (Q + t(Q))/2
res <- run_experiments(Q,
                       ns=ns,
                       sim=sim,
                       mc_cores=mc_cores,
                       nlambda = nlambda,
                       lambda.min.ratio = lambda.min.ratio,
                       alpha.grid=alpha.grid)

out <- summarize_plot_results(res)

# view the summary table
out$table

# display the RMSE plot
print(out$plots$rmse_grid)

# display the false‐discovery plot
print(out$plots$rate_grid)
