### Kullback-Leibler Divergence in Simulation Study
library(MASS)
library(kldest)

load('Post1D.RData')
load('AEArun1D.RData')
AEA <- result$beta[401:60000] #delete burn-in

kl_divergence_mcmc <- function(chain_p, chain_q, n_grid = 2048, eps = 1e-12) {
  # Kernel density estimates
  dens_p <- density(chain_p, n = n_grid)
  dens_q <- density(chain_q, n = n_grid,
                    from = min(dens_p$x),
                    to   = max(dens_p$x))
  
  # Interpolate q onto p's grid
  q_interp <- approx(dens_q$x, dens_q$y, xout = dens_p$x)$y
  
  # Avoid log(0)
  p <- dens_p$y + eps
  q <- q_interp + eps
  
  # Grid spacing
  dx <- diff(dens_p$x)[1]
  
  # KL divergence
  kl <- sum(p * log(p / q)) * dx
  return(kl)
}


KL_matrix <- matrix(0,4,5)
for(i in 1:5){
  KL_matrix[1,i] <- kl_divergence_mcmc(AEA, PostAll[[i]])
  KL_matrix[2,i] <- kl_divergence_mcmc(AEA, PostAll[[i+5]])
  KL_matrix[3,i] <- kl_divergence_mcmc(AEA, PostAll[[i+10]])
  KL_matrix[4,i] <- kl_divergence_mcmc(AEA, PostAll[[i+15]])
}

PostMean <-lapply(PostAll, mean)
Post_Mean <- matrix(0,4,5)
Post_Mode <- matrix(0,4,5)
densAll <- lapply(PostAll, density)
for (i in 1:5){
  Post_Mean[1,i] <- PostMean[[i]]
  Post_Mean[2,i] <- PostMean[[i+5]]
  Post_Mean[3,i] <- PostMean[[i+10]]
  Post_Mean[4,i] <- PostMean[[i+15]]
  Post_Mode[1,i] <- densAll[[i]]$x[which.max(densAll[[i]]$y)]
  Post_Mode[2,i] <- densAll[[i+5]]$x[which.max(densAll[[i+5]]$y)]
  Post_Mode[3,i] <- densAll[[i+10]]$x[which.max(densAll[[i+10]]$y)]
  Post_Mode[4,i] <- densAll[[i+15]]$x[which.max(densAll[[i+15]]$y)]
}

allPap <- matrix(0,4,10)
for(i in 1:5){
  allPap[,2*i-1] <- Post_Mean[,i]
  allPap[,2*i] <- KL_matrix[,i]
}

rownames(KL_matrix)<- c('EqLin', 'EqHer', 'ExpLin', 'ExpHer')
colnames(KL_matrix)<- c('5','10','20','40', '60')
rownames(Post_Mean)<- c('EqLin', 'EqHer', 'ExpLin', 'ExpHer')
colnames(Post_Mean)<- c('5','10','20','40', '60')
rownames(Post_Mode)<- c('EqLin', 'EqHer', 'ExpLin', 'ExpHer')
colnames(Post_Mode)<- c('5','10','20','40', '60')

