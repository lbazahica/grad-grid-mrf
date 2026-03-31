## AEA for Menteith data

library(bayess)
data("Menteith")
iter <- 1600
burn <- iter/4 + 1
n <- prod(dim(Menteith))
k <- 6
image(as.matrix(Menteith),asp=1,xaxt='n',yaxt='n',col=gray(0:255/255))


priors <- list()
priors$k <- k
priors$mu <- seq(min(Menteith),max(Menteith),by=diff(range(Menteith))/5)
priors$mu.sd <- rep(diff(range(Menteith))/20,k)
priors$sigma <- rep(diff(range(Menteith))/20,k)
priors$sigma.nu <- rep(3, k)
priors$beta <- c(0,2)

samp.mu <- rnorm(k*iter, 0, priors$mu.sd[1])
samp.mu <- matrix(samp.mu + rep(priors$mu, each=iter), ncol=k, nrow=iter)
colMeans(samp.mu)
apply(samp.mu,2,sd)
old.par <- par(mfrow=c(3,2))
for (i in 1:k) {
  hist(samp.mu[,i], breaks=20, main=bquote(mu[.(i)]), xlab=bquote(mu[.(i)]), freq=F)
  curve(dnorm(x,mean=priors$mu[i], sd=priors$mu.sd[i]), col=4, lty=2, lwd=2, add=T)
}
par(old.par)
hist(samp.mu, breaks=70, freq=FALSE, xlab=expression(mu), main="", xaxt='n')
for (i in 1:k) {
  curve(dnorm(x,mean=priors$mu[i], sd=priors$mu.sd[i])/k, col=i, lty=i, lwd=2, add=T)
}
axis(1, at=priors$mu)
title(main="Priors for means")


SS0 <- priors$sigma.nu*priors$sigma^2
s0 <- 1/sqrt(rgamma(iter,priors$sigma.nu/2,SS0/2))
hist(s0,breaks=20,freq=F,xlab=expression(sigma),
     main="Prior for standard deviation")
lines(density(s0), lwd=3, lty=2, col=4)


range(Menteith)
y <- as.numeric(unlist(Menteith))
hist(y, freq=F, xlab="", main="Pixel Values", breaks=20, xaxt='n')
axis(1, at=priors$mu)
lines(density(y), lwd=3, lty=2, col=4)


library(bayesImageS)
bcrit <- log(1 + sqrt(k))
mask <- matrix(1, nrow=sqrt(n), ncol=sqrt(n))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
E0 <- maxS/k
V0 <- maxS*(1/k)*(1 - 1/k)
it <- 6e4
burn <- 401
mh <- list(algorithm="ex", bandwidth=0.02, auxiliary=400)
set.seed(1234)
tm <- Sys.time()
result <- mcmcPotts(y,neigh,block,priors,mh,it,burn)
result$time <- Sys.time() - tm
save(result, file="result.rda")
print(paste(format(result$time), "for",length(result$beta),"MCMC iterations."))
print(paste("MH acceptance rate",format(result$accept/it),"for algorithm",mh$algorithm))