# Load required package
library(dplyr)

# Define parameters of population model II
parasample_pop <- expand.grid(
  beta_0 = list(-1.5, -0.43),                       # Latent class intercepts 
  beta_x1x2 = list(c(-0.50, 0.75)),                       # Latent class slopes (covariates x1 and x2)
  alpha_0 = list(c(1.5, -1.5), c(0.5, -0.5)),        # Latent intercept factor 
  gamma_0 = list(c(0.5, -0.5)),                      # Latent intercept factor i1
  alpha_1 = list(c(-0.2, 0.2)),                      # Latent slope factor
  theta = list(rep(1, 18)),                       # Variance of items for 6 time points
  sigma = list(c(1, -0.15, 1)),                   # Variance-covariance for class 1
  sigma = list(c(1, -0.15, 1))                    # Variance-covariance for class 2
) %>% 
  `rownames<-`(c("EH", "UH", "EL", "UL"))           # Define row names using %>%