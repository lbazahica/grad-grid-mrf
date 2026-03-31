## AEA run for simulation study

library(bayesImageS)

set.seed(1) 
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
k <- 6
beta.true <- 1.274
res.sw <- swNoData(beta.true, k, neigh, blocks, niter=200)
z <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))

## Prior configuration, need to test this more strategically
priors <- list()
priors$k <- k
priors$mu <- 2*c(-5,-3,-1,1,3,5) #c(-4, -2, 0, 2, 4, 6)
#priors$mu <- c(-4, -2, 0, 2, 4) #for k=5
priors$mu.sd <- rep(0.25,k)
priors$sigma <- rep(0.125,k)
priors$sigma.nu <- rep(3, k)
priors$beta <- c(0,2.5)

## Check values for m0 and s0
m0 <- sort(rnorm(priors$k,priors$mu,priors$mu.sd))
SS0 <- priors$sigma.nu*priors$sigma^2
s0 <- 1/sqrt(rgamma(priors$k,priors$sigma.nu/2,SS0/2))
l <- as.vector(z)
y <- m0[l] + rnorm(nrow(neigh),0,s0[l])
#save(y, file="SimulationData.RData")

#AEA
iter <- 6e4 #needs to be 6e4, runs approx. 10 hours
burn <- 401 
#mh <- list(algorithm="ex",bandwidth=0.2,auxiliary=500)#, aux_swap=TRUE) 
mh <- list(algorithm="ex", bandwidth=0.02, auxiliary=400)
tm <- system.time(result <- mcmcPotts(y,neigh,blocks,priors,mh,iter,burn))
result$time <- tm
print(paste(format(result$time[3]), "seconds for",length(result$beta),"MCMC iterations."))
print(paste("MH acceptance rate",format(result$accept/iter),"for algorithm",mh$algorithm))


#Visualisation
plot.ts(result$beta[1:(2*burn)], ylab=expression(beta), xlab="Iteration")
abline(v=burn,col=2,lty=3)
hist(result$beta[burn:iter], breaks=20, xlab=expression(beta),
     main="MCMC samples for inverse temperature", freq=F)
lines(density(result$beta[burn:iter], bw=0.001), lwd=3, lty=2, col=4)

#save(result, file = "AEArun1D.RData")