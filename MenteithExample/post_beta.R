## 1. precomputation step ----
library(bayesImageS)
data(Menteith)
mask <- matrix(1,nrow(Menteith),ncol(Menteith))
n <- prod(dim(mask))
k <- 6
bcrit <- log(1 + sqrt(k))
beta <- seq(0,2.5,by=0.1)
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
E0 <- maxS/k
V0 <- maxS*(1/k)*(1 - 1/k)

library(doParallel)
cores <- min(detectCores(), length(beta), 4)

# Gradient Method ----

library(mcmcse)
min_step <- 0.1
iter <- 800
burn <- iter/4 + 1

beta_g <- round(bcrit, digits=2)
res <- swNoData(beta_g,k,neigh,block,iter)
matu_g <- res$sum
beta_g_chain <- res$sum[burn:iter]
beta_g_ess <- ceiling((iter-burn+1)/ess(beta_g_chain))
beta_g_esschain <- beta_g_chain[seq(1, length(beta_g_chain), beta_g_ess)]
beta_g_var <- var(beta_g_esschain)
min_step <- exp(-4.25)
beta_g3 <- c(beta_g-min_step, beta_g, beta_g+min_step)
min_step <- (1+beta_g_var)^(-1/8)

source("calc_beta.R")

# New Gradient 1 ----
matu_g3 <- matu_g
ExpressionVect <- c(substitute(dir_low(beta_g3[1], 2)), 
                    substitute(dir_high(beta_g3[length(beta_g3)], 2)))

tm3 <- system.time(results <- mclapply(ExpressionVect, eval, mc.cores= 2))
print(tm3)
beta_g3 <- c(results[[1]]$beta_low, beta_g3[2], results[[2]]$beta_high)
matu_g3 <- cbind(results[[1]]$matu_low, matu_g3, results[[2]]$matu_high)

getSubSamples <- function(beta, matu){
  iter <- nrow(matu)
  burn <- iter/4+1
  subSample <- ess(matu[burn:iter,])
  simEss <-  matrix(nrow=3, ncol=length(beta))
  rownames(simEss) <- c("beta","mu", "sd")
  for(i in 1:length(beta)){
    simEss[1,i] <- beta[i]
    subS <- sample(matu[burn:iter,i], round(subSample[i]))
    simEss[2,i] <- mean(subS)
    simEss[3,i] <- sd(subS)
  }
  return(simEss)
}

g3Ess <- getSubSamples(beta_g3, matu_g3)
## Linear Interpolation log
lrcst3=approxfun(beta_g3,g3Ess[2,])
#hermite log interpolation
library(signal)
m <- 10*length(beta)
range <- seq(min(beta), max(beta), by = 1e-5)
hermPolg3<- pchip(beta_g3,g3Ess[2,], range)

beta_g <- list(beta=beta_g3, matu=matu_g3, linfnc=lrcst3, hermPol=hermPolg3)
save(beta_g, file = "betaExp1.RData")