### Exchange algorithm with the hidden Potts model

# Load source
library(bayess)
library(bayesImageS)

data(Menteith)
load('MenteithBeta.RData') #renamed Beta from Simulation Study
linfnc <-  MenteithBeta$betaEqLin10

# Setup from pixel data
set.seed(1234)
mask <- matrix(1,nrow(Menteith),ncol(Menteith))
n <- prod(dim(mask))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
k <- 6
bcrit <- log(1 + sqrt(k))
y <- as.matrix(Menteith)

# Priors
priors <- list()
priors$k <- k
priors$mu <- seq(min(Menteith),max(Menteith),by=diff(range(Menteith))/5)
priors$mu.sd <- rep(diff(range(Menteith))/20,k)
priors$sigma <- rep(diff(range(Menteith))/20,k)
priors$sigma.nu <- rep(3, k)
priors$beta <- c(0,2)

### exchange algorithm
## burn-in
nMC <- 60000 # 2201
nburn <- 401
beta <- rep(0, nMC)
beta[1] <- 1
acc <- 0
sumZs <- numeric(nMC)
muEst <- sdEst <- matrix(nrow=nMC, ncol=priors$k)

# initialize the parameters by sampling from the prior
prSS <- priors$sigma.nu * priors$sigma^2
muEst[1,] <- rnorm(priors$k, priors$mu, priors$mu.sd)
sdEst[1,] <- 1/sqrt(rgamma(priors$k, priors$sigma.nu/2, prSS/2))

# initialize the pixel labels 
n <- length(y)
z <- matrix(0, nrow=n+1, ncol=priors$k)
randZ <- sample.int(priors$k, n, replace = TRUE)
for (i in 1:n) {
  z[i,(randZ[i])] <- 1
}

set.seed(12)

t1 <- proc.time()
for (i in 2:nMC){

  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, block, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Simulate
  summary_stats_star <- linfnc(beta_star)
  
  log_MH_ratio <- beta_star * (summary_stats - summary_stats_star) +
    beta[i-1] * (summary_stats_star - summary_stats)
  log_accept_prob <- log(runif(1, 0, 1))
  if (log_MH_ratio > log_accept_prob){
    beta[i] <- beta_star
    if (i >= nburn) {acc <- acc + 1}
  }
  else{
    beta[i] <- beta[i-1]
  }
  #print(paste(i, beta[i], sumZs[i], acc))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

pathAlg <- list(beta, muEst, sdEst, sumZs, z, acc, time_total, time_burnin)
save(pathAlg, file = 'hidden_potts_lin.rda')
