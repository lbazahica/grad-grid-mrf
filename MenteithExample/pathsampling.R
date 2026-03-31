## Path sampling/thermodynamic integration for Menteith example

library(bayesImageS)

data("Menteith")
set.seed(1234)
mask <- matrix(1,100,100)
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
k <- 6
beta <- seq(0,2.5,length.out=10)
iter <- 1600 
burn <- iter/4 + 1
n <- prod(dim(mask))
y <- as.numeric(unlist(Menteith))

priors <- list()
priors$k <- k
priors$mu <- seq(min(Menteith),max(Menteith),by=diff(range(Menteith))/5)
priors$mu.sd <- rep(diff(range(Menteith))/20,k)
priors$sigma <- rep(diff(range(Menteith))/20,k)
priors$sigma.nu <- rep(3, k)
priors$beta <- c(0,2)

library(doParallel)
cores <- min(detectCores(), length(beta), 2)
print(paste("Parallel computation using",cores,"CPU cores:",
            iter,"iterations for",length(beta),"values of beta."))
cl <- makeCluster(cores)
print(cl)
clusterSetRNGStream(cl)
registerDoParallel(cl)

tm <- system.time(matu <- foreach(i=1:length(beta),
                                  .packages=c("bayesImageS"), .combine='cbind') %dopar% {
                                    res <- swNoData(beta[i],k,neigh,block,iter)
                                    res$sum
                                  })
stopCluster(cl)
#data("simSW", package = "bayesImageS")
simSW <- list(matu = matu, tm = tm)
print(simSW$tm)

## thermodynamic integration (path sampling) ----
pathMx <- matrix(nrow=2, ncol=length(beta))
rownames(pathMx) <- c("beta","mu")
pathMx[1,] <- beta
pathMx[2,] <- colMeans(simSW$matu)
iter <- 60000
burn <- 401

alg <- list()
alg$path <- list(algorithm="path",bandwidth=0.02,path=pathMx)
mh <- alg$path
tm <- system.time(result <- mcmcPotts(y,neigh,block,priors,mh,iter,burn))
print(tm["elapsed"])
#print(res.ex$elapsed)
result$elapsed <- tm["elapsed"]
result$system <- tm["sys.self"]
result$precomp <- simSW$tm
pathData <- result
save('pathData', file='hidden_potts_path.rda')
