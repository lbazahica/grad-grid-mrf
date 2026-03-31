# Simulation Study in 1D
library(bayesImageS)
library(bayess)
library(mcmcse)

set.seed(1234)

## 1. precomputation step ----
k <- 6
bcrit <- log(1 + sqrt(k))
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
E0 <- maxS/k

iter <- 1600
burn <- iter/4 + 1
n <- prod(dim(mask))

# Gradient Methods ----
faktor <- c(-1,-4.4,-8.5,-12.05,-16.76)

beta_g <- round(bcrit, digits=2)
res <- swNoData(beta_g,k,neigh,block,iter)
matu_g <- res$sum
beta_g_chain <- res$sum[burn:iter]
beta_g_ess <- ceiling((iter-burn+1)/ess(beta_g_chain))
beta_g_esschain <- beta_g_chain[seq(1, length(beta_g_chain), beta_g_ess)]
beta_g_var <- var(beta_g_esschain)

source('calc_beta.R')
library(parallel)
expgrids <- 0
for(j in 1:length(faktor)) {
  min_step <- exp(faktor[j])
  beta_g1 <- c(beta_g - min_step, beta_g, beta_g + min_step)
  
  # Gradient with exp-term ----
  matu_g1 <- matu_g
  ExpressionVect <- c(substitute(dir_low(beta_g1[1], step = 2, faktor = faktor[j])),
                      substitute(dir_high(beta_g1[length(beta_g1)], step = 2, faktor = faktor[j])))
  
  tm1 <- system.time(results <- mclapply(ExpressionVect, eval, mc.cores = 2))
  print(tm1)
  beta_g1 <- c(results[[1]]$beta_low, beta_g1[2], results[[2]]$beta_high)
  matu_g1 <- cbind(results[[1]]$matu_low, matu_g1, results[[2]]$matu_high)
  expgrids <- c(expgrids, list(beta_g1, matu_g1, tm1))
}
print(length(beta_g1))
var <- c("beta", "matu", "tm")
expgrids[[1]] <- NULL 
names(expgrids) <- paste0(var, "_", rep(c(5,10,20,40,60),each=3))
# save(expgrids, file = "expgridset1D.RData")
