### Summarize the results of the example with the hidden Potts model

# Load source
library(ggplot2)

## Data
load('hidden_potts_path.rda')
load('hidden_potts_lin.rda')
load('hidden_potts_exp.rda')
load("result.rda") #AEA data


## Path sampling
beta_PA <- pathData$beta[401:length(pathData[[1]])]
# Posterior mean
post_mean_PA <- mean(beta_PA)
# Posterior sd
post_sd_PA <- sd(beta_PA)
# Time
time_PA <- as.numeric(pathData$elapsed)/3600
# ESS/time
ESS_time_PA <- coda::effectiveSize(beta_PA) / time_PA

## Equidistant with linear Interpolation
beta_DA <- pathAlg[[1]][401:length(pathAlg[[1]])]
# Posterior mean
post_mean_DA <- mean(beta_DA)
# Posterior sd
post_sd_DA <- sd(beta_DA)
# Time
time_DA <- (pathAlg[[7]][3] - pathAlg[[8]][3]) / 3600
# ESS/time
ESS_time_DA <- coda::effectiveSize(beta_DA) / time_DA

## Gradient-based with hermite interpolation
beta_EP <- exchange[[1]][401:length(exchange[[1]])]
# Posterior mean
post_mean_EP <- mean(beta_EP)
# Posterior sd
post_sd_EP <- sd(beta_EP)
# Time
time_EP <- (exchange[[7]][3] - exchange[[8]][3]) / 3600
# ESS/time
ESS_time_EP <- coda::effectiveSize(beta_EP) / time_EP

## Approximate Exchange Algorithm
beta_EX <- result$beta[401:length(result$beta)]
# Posterior mean
post_mean_EX <- mean(beta_EX)
# Posterior sd
post_sd_EX <- sd(beta_EX)
# Time
time_EX <- as.numeric(result$time)
# ESS/time
ESS_time_EX <- coda::effectiveSize(beta_EX) / time_EX

## Plot
df_PA <- data.frame(beta = beta_PA)
df_DA <- data.frame(beta = beta_DA)
df_EP <- data.frame(beta = beta_EP)
df_EX <- data.frame(beta = beta_EX)

plot2 <- ggplot() + 
  stat_density(data = df_PA, aes(beta, color = "TI"), adjust = 2, geom = "line", linetype = "dashed", linewidth = 0.8) +
  stat_density(data = df_EP, aes(beta, color = "GH"), adjust = 2, geom = "line", linetype = "dashed", linewidth = 0.8) +
  stat_density(data = df_DA, aes(beta, color = "EL"), adjust = 2, geom = "line", linetype = "dashed", linewidth = 0.8) +
  stat_density(data = df_EX, aes(beta, color = "AEA"), adjust = 2, geom = "line", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("TI" = "#e41a1c", "EL" = "#377eb8", "GH" = "#4daf4a", "AEA" = "#984ea3")) +
  xlim(c(1.23, 1.38)) +
  theme_bw() +
  labs(x = expression(beta^"(1)"), y = "density", color = "Approach") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )

plot2

#ggsave(filename="hidden_potts_posterior_density_plot.png", plot=plot2, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

