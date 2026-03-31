### Create 1-dimensional testset for simulation study

library(bayesImageS)
library(bayess)

set.seed(1234)
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
k <- 6
block <- getBlocks(mask, 2)
n <- prod(dim(mask))

## 1. precomputation step ----
beta_mid <- sort(runif(1000, 0, 2.5))
iter_mid <- 1600
burn_mid <- iter_mid/4+1
library(doParallel)
cores <- min(detectCores(), length(beta_mid), 4)
print(paste("Parallel computation using",cores,"CPU cores:",
            iter_mid,"iterations for",length(beta_mid),"values of beta."))
cl <- makeCluster(cores)
print(cl)
clusterSetRNGStream(cl)
registerDoParallel(cl)

tm_mid <- system.time(matu_mid <- foreach(i=1:length(beta_mid),
                                          .packages=c("bayesImageS"), .combine='cbind') %dopar% {
                                            res <- swNoData(beta_mid[i],k,neigh,block,iter_mid)
                                            res$sum
                                          })

stopCluster(cl)
testset <- list(beta_mid, matu_mid, tm_mid)
save(testset, file = "testset.RData")
