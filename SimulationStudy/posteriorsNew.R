## Posteriors for simulation study data

library(bayesImageS)
library(bayess)

#load interpolated grid data (from file created in RMSEall1D.R)
load('betaInt.Rdata')

#calculate posteriors
set.seed(10) 
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
k <- 6
beta.true <- 1.274

n <- prod(dim(mask))
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
bcrit <- log(1 + sqrt(k))

# Priors
priors <- list()
priors$k <- k
priors$mu <- 2*c(-5,-3,-1,1,3,5)
priors$mu.sd <- rep(0.25,k)
priors$sigma <- rep(0.125,k)
priors$sigma.nu <- rep(3, k)
priors$beta <- c(0,2.5)

load('SimulationData.RData')

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

## MCMC calculations ----
betaPostEq <- list(1,2,3,4,5)
betaPostExp <- list(1,2,3,4,5)
betaPostEqLin <- list(1,2,3,4,5)
betaPostExpLin <- list(1,2,3,4,5)
timesEqHer <- rep(0,5)
timesExpHer <- rep(0,5)
timesEqLin <- rep(0,5)
timesExpLin <- rep(0,5)

for(j in 1:5){
  ## equidistant grid , hermite interpolation ----
hermPol <- betaEqHer[[j+1]]
t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, blocks, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Interpolate
  hilfe <- round(beta_star, 5)*1e5
  summary_stats_star <- hermPol[hilfe]
  
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
betaPostEq[[j]] <- beta
timesEqHer[j] <- time_total[3]

## equidistant grid , linear interpolation ----
linfnc <- betaEqLin[[j+1]]
t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, blocks, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Interpolate
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
betaPostEqLin[[j]] <- beta
timesEqLin[j] <- time_total[3]

## gradient-based grid , hermite interpolation ----
hermPol <- betaExpHer[[j+1]]
t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, blocks, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Interpolate
  hilfe <- round(beta_star, 5)*1e5
  summary_stats_star <- hermPol[hilfe]
  
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
betaPostExp[[j]] <- beta
timesExpHer[j] <- time_total[3]

## gradient-based grid , linear interpolation ----
linfnc <- betaExpLin[[j+1]]
t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, blocks, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Interpolate
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
betaPostExpLin[[j]] <- beta
timesExpLin[j] <- time_total[3]
}

### Gathering posterior data ----
df_1 <- data.frame(beta = betaPostEq[[1]][nburn:nMC])
df_2 <- data.frame(beta = betaPostEq[[2]][nburn:nMC])
df_3 <- data.frame(beta = betaPostEq[[3]][nburn:nMC])
df_4 <- data.frame(beta = betaPostEq[[4]][nburn:nMC])
df_5 <- data.frame(beta = betaPostEq[[5]][nburn:nMC])
df_6 <- data.frame(beta = betaPostExp[[1]][nburn:nMC])
df_7 <- data.frame(beta = betaPostExp[[2]][nburn:nMC])
df_8 <- data.frame(beta = betaPostExp[[3]][nburn:nMC])
df_9 <- data.frame(beta = betaPostExp[[4]][nburn:nMC])
df_10 <- data.frame(beta = betaPostExp[[5]][nburn:nMC])

dfL_1 <- data.frame(beta = betaPostEqLin[[1]][nburn:nMC])
dfL_2 <- data.frame(beta = betaPostEqLin[[2]][nburn:nMC])
dfL_3 <- data.frame(beta = betaPostEqLin[[3]][nburn:nMC])
dfL_4 <- data.frame(beta = betaPostEqLin[[4]][nburn:nMC])
dfL_5 <- data.frame(beta = betaPostEqLin[[5]][nburn:nMC])
dfL_6 <- data.frame(beta = betaPostExpLin[[1]][nburn:nMC])
dfL_7 <- data.frame(beta = betaPostExpLin[[2]][nburn:nMC])
dfL_8 <- data.frame(beta = betaPostExpLin[[3]][nburn:nMC])
dfL_9 <- data.frame(beta = betaPostExpLin[[4]][nburn:nMC])
dfL_10 <- data.frame(beta = betaPostExpLin[[5]][nburn:nMC])

PostAll <- list(EqLin5 = dfL_1$beta, EqLin10=dfL_2$beta, EqLin20 = dfL_3$beta, EqLin40=dfL_4$beta, EqLin80=dfL_5$beta,
                EqHer5 = df_1$beta, EqHer10=df_2$beta, EqHer20 = df_3$beta, EqHer40=df_4$beta, EqHer80=df_5$beta,
                ExpLin5 = dfL_6$beta, ExpLin10=dfL_7$beta, ExpLin20 = dfL_8$beta, ExpLin40=dfL_9$beta, ExpLin60=dfL_10$beta,
                ExpHer5 = df_6$beta, ExpHer10=df_7$beta, ExpHer20 = df_8$beta, ExpHer40=df_9$beta, ExpHer60=df_10$beta,
                timesEqLin = timesEqLin, timesEqHer= timesEqHer, timesExpLin=timesExpLin, timesExpHer=timesExpHer)

save('PostAll', file='Post1D.RData')
