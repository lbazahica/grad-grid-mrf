calc_beta <- function(beta, step = 0, min_step = 0.1){
  res <-swNoData(beta,k,neigh,block,iter)
  chain <- res$sum[burn:iter]
  ess <- ceiling((iter-burn+1)/ess(chain))
  esschain <- chain[seq(1, length(chain), ess)]
  if(step == 0){
    stepsize <- 0.5*sqrt(V0)/sd(esschain)#/3
  }else if(step == 1){
    stepsize <- min_step+(1-min_step)*(1-log(var(esschain))/log(beta_g_var))
  }else if(step == 2){
    stepsize <- exp(-4.37*var(esschain)/beta_g_var)
  }else if(step == 3){
    stepsize <- (1+var(esschain))^(-1/8)
  }else if(step == 4){
    stepsize <- (1+var(esschain))^(-1/12)*exp(-2*var(esschain)/beta_g_var)
  }else if(step == 5){
    stepsize <- (1+beta_g_var+var(esschain))^(-21/220)
  }else{
    stop("unsuitable step choice")
  }
  return(list(stepsize = stepsize, chain = res$sum))
}

dir_low <- function(beta, step = 0, min_step = 0.1){
  matu_low <- matrix(, nrow = iter, ncol = 0)
  while(beta[1]!=0){
    beta_low <- calc_beta(beta[1], step, min_step)
    matu_low <- cbind(beta_low$chain, matu_low)
    beta_next <- beta[1]-beta_low$stepsize
    if(beta_next>0){
      beta <- c(beta_next, beta)
    }else{
      beta <- c(0, beta)
    }
  }
  beta_low <- calc_beta(beta[1], step, min_step)
  matu_low <- cbind(beta_low$chain, matu_low)
  return(list(beta_low = beta, matu_low = matu_low))
}

dir_high <- function(beta, step = 0, min_step = 0.1){
  matu_high <- matrix(, nrow = iter, ncol = 0)
  while(beta[length(beta)]!=2.5){
    beta_high <- calc_beta(beta[length(beta)], step, min_step)
    matu_high <- cbind(matu_high, beta_high$chain)
    beta_next <- beta[length(beta)]+beta_high$stepsize
    if(beta_next<2.5){
      beta <- c(beta, beta_next)
    }else{
      beta <- c(beta, 2.5)
    }
  }
  beta_high <- calc_beta(beta[length(beta)], step, min_step)
  matu_high <- cbind(matu_high, beta_high$chain)
  return(list(beta_high = beta, matu_high = matu_high))
}