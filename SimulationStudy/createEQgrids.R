# Simulation Study in 1D
# create equidistant grids
library(bayesImageS)
library(bayess)
library(mcmcse)
N <- c(5,10,20,40,60)

set.seed(1234)
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
k <- 6

iter <- 1600
burn <- iter/4 + 1
n <- prod(dim(mask))

## 1. precomputation step ----
bcrit <- log(1 + sqrt(k))
mask <- matrix(1, nrow=sqrt(n), ncol=sqrt(n))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
E0 <- maxS/k
V0 <- maxS*(1/k)*(1 - 1/k)

library(doParallel)
cores <- min(detectCores(), length(beta), 4)
print(paste("Parallel computation using",cores,"CPU cores:",
            iter,"iterations for",length(N),"sequences of beta values."))
cl <- makeCluster(cores)
print(cl)
clusterSetRNGStream(cl)
registerDoParallel(cl)

gridset <- list()
tm <- 0
for(j in 1:length(N)){
  beta <- seq(0,2.5,length=N[j])
  tm <- system.time(matu <- foreach(i=1:length(beta),
                                  .packages=c("bayesImageS"), .combine='cbind') %dopar% {
                                    res <- swNoData(beta[i],k,neigh,block,iter)
                                    res$sum
                                  })
  print(paste('step', j, 'done'))
  gridset <- c(gridset, list(beta, matu, tm))
}
var <- c("beta", "matu", "tm")
names(gridset) <- paste0(var, "_", rep(N[1:5],each=3))
save(gridset, file = "eqgridset1D.RData")
