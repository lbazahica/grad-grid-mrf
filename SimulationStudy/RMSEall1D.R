# Simulation Study in 1D -- all RMSEs
library(mcmcse)
library(signal)
load('testset.RData') #beta in [0,2.5]
load('eqgridset1D.RData')
load('expgridset1D.RData')

set.seed(36)
iter <- 1600 
burn <- iter/4 + 1

## Subsample from Testset ----
iter_mid <- nrow(testset[[2]])
burn_mid <- iter_mid/4+1
beta_mid <- testset[[1]]
matu_mid <- testset[[2]]
simSW_mid <- list(matu = matu_mid)

subSample <- pmin(ess(simSW_mid$matu[burn_mid:iter_mid,]), iter_mid-burn_mid+1)
simEss_mid <-  matrix(nrow=3, ncol=length(beta_mid))
rownames(simEss_mid) <- c("beta","mu", "sd")
for(i in 1:length(beta_mid)){
  simEss_mid[1,i] <- beta_mid[i]
  subS <- sample(simSW_mid$matu[burn_mid:iter_mid,i], round(subSample[i]))
  simEss_mid[2,i] <- mean(subS)
  simEss_mid[3,i] <- sd(subS)
}


## Equidistant set ----

N <- c(5,10,20,40,60)
RMSE <- matrix(0, nrow=length(N), ncol = 2)
colnames(RMSE) <- c("Linear", "Hermite")
betaEqLin <- list(N=N)
betaEqHer <- list(N=N)

for(j in 1:length(N)) {
  beta <- seq(0, 2, length = N[j])
  var <- c("beta", "matu")
  var2 <- paste0(var, "_", N[j])
  beta <- gridset[[var2[1]]]
  matu <- gridset[[var2[2]]]
  
  simSW <- list(matu = matu)#, tm = tm)
  #print(simSW$tm)
  
  subSample <- ess(simSW$matu[burn:iter, ])
  simEss <-  matrix(nrow = 3, ncol = length(beta))
  rownames(simEss) <- c("beta", "mu", "sd")
  for (i in 1:length(beta)) {
    simEss[1, i] <- beta[i]
    subS <- sample(simSW$matu[burn:iter, i], round(subSample[i]))
    simEss[2, i] <- mean(subS)
    simEss[3, i] <- sd(subS)
  }
  
  
  ## Equidistant grid 
  
  ## Linear Interpolation equidistant
  lrcst = approxfun(beta, simEss[2, ])
  betaEqLin <- append(betaEqLin, lrcst)
  names(betaEqLin)[j+1] <- var2[1]
  
  #Residuals linear
  res_L <- lrcst(beta_mid) - simEss_mid[2, ]
  
  ## Hermite Interpolation
  m <- 10 * length(beta)
  range <- seq(min(beta), max(beta), by = 1e-5)
  hermPol <- pchip(beta, simEss[2, ], range)
  betaEqHer <- append(betaEqHer, list(hermPol))
  names(betaEqHer)[j+1] <- var2[1]
  
  #Residuals hermite
  hilfe <- c(round(beta_mid, 5) * 1e5)
  compH <- hermPol[hilfe]
  res_H <- compH - simEss_mid[2, ]
  
  RMSE[j, ] <- c(sqrt(sum(abs(res_L)^2)/length(abs(res_L))), sqrt(sum(abs(res_H)^2) /
                                                                    length(abs(res_H))))
}


## Gradient-based set ----

N2 <- c(5,10,20,40,60)
RMSEEXP <- matrix(0, nrow=length(N2), ncol = 2)
colnames(RMSEEXP) <- c("Linear", "Hermite")
expgrids$beta_60[16] <- expgrids$beta_60[16]+0.0000000000001
expgrids$beta_60[27] <- expgrids$beta_60[27]+0.0000000000001
expgrids$beta_60[36] <- expgrids$beta_60[36]-0.0000000000002
expgrids$beta_60[37] <- expgrids$beta_60[37]-0.0000000000001
expgrids$beta_60[39] <- expgrids$beta_60[39]+0.0000000000001

betaExpLin <- list(N=N2)
betaExpHer <- list(N=N2)

for(j in 1:length(N2)) {
  beta <- seq(0, 2, length = N2[j])
  var <- c("beta", "matu")
  var2 <- paste0(var, "_", N2[j])
  beta <- expgrids[[var2[1]]]
  matu <- expgrids[[var2[2]]]
  
  simSW <- list(matu = matu)#, tm = tm)
  #print(simSW$tm)
  
  subSample <- ess(simSW$matu[burn:iter, ])
  simEss <-  matrix(nrow = 3, ncol = length(beta))
  rownames(simEss) <- c("beta", "mu", "sd")
  for (i in 1:length(beta)) {
    simEss[1, i] <- beta[i]
    subS <- sample(simSW$matu[burn:iter, i], min(length(simSW$matu[burn:iter, i]), round(subSample[i])))
    simEss[2, i] <- mean(subS)
    simEss[3, i] <- sd(subS)
  }
  
  
  ## Exponential grid 
  
  ## Linear Interpolation equidistant
  lrcst = approxfun(beta, simEss[2, ])
  betaExpLin <- append(betaExpLin, lrcst)
  names(betaExpLin)[j+1] <- var2[1]
  
  #Residuals linear
  res_L <- lrcst(beta_mid) - simEss_mid[2, ]
  
  ## Hermite Interpolation
  m <- 10 * length(beta)
  range <- seq(min(beta), max(beta), by = 1e-5)
  hermPol <- pchip(beta, simEss[2, ], range)
  betaExpHer <- append(betaExpHer, list(hermPol))
  names(betaExpHer)[j+1] <- var2[1]
  
  #Residuals hermite
  hilfe <- c(round(beta_mid, 5) * 1e5)
  compH <- hermPol[hilfe]
  res_H <- compH - simEss_mid[2, ]
  
  RMSEEXP[j, ] <- c(sqrt(sum(abs(res_L)^2)/length(abs(res_L))), sqrt(sum(abs(res_H)^2) /
                                                                       length(abs(res_H))))
}

## Plots ----

y <- c(min(RMSE, RMSEEXP)-10, max(RMSE, RMSEEXP)+10)
plot(x=N, y=RMSE[,1],xlab='number of grid points',ylab=expression(RMSE),
     type='b',pch=15, lty=2, log='xy', ylim=y, cex.lab = 1.2)
points(N, RMSE[,2], col="blue", type='b', lty=2, pch=16)
points(N2, RMSEEXP[,1], col="red", type='b',lty=2, pch=17)
points(N2, RMSEEXP[,2], col="darkgreen", type='b',lty=2, pch=18)
#title(main=paste("RMSE for different grid sizes"))
legend("topright", legend=c("EL","EH","GL","GH"), 
       pch=c(15,16,17,18),col=c("black", "blue","red","darkgreen"), lty=2)

save(list=c('betaEqHer', 'betaEqLin', 'betaExpHer', 'betaExpLin'), file='betaInt.Rdata')
