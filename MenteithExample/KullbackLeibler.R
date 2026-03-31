### Kullback-Leibler Divergence
library(MASS)

# load data
load("hidden_potts_path.rda")
path <- pathData$beta[401:length(pathData$beta)]

load("hidden_potts_lin.rda")
equidistant <- pathAlg[[1]][401:length(pathAlg[[1]])]

load("hidden_potts_exp.rda")
exponential <- exchange[[1]][401:length(exchange[[1]])]

load("/Users/z142094/Documents/Potts/Final/MenteithExample/result.rda")
exchange <- result$beta[401:length(result$beta)]

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

kl_ti_ex <- kl_divergence_mcmc(path,exchange)
kl_ex_ex <- kl_divergence_mcmc(exponential,exchange)
kl_eq_ex <- kl_divergence_mcmc(equidistant,exchange)

