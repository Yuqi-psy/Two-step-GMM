
set.seed(1119)
Ninds <- 3 # set the number of indicators
Nobs <- 3 # repeated measures
Nclass <- 2 # number of observed variables and number of latent classes
Nsubs <- c(5000,5000) # number of subjects in each latent class (equal)
Npops <- 10000 # number of the whole populations
Cov_gf <- c(0,0.3) # the correlation of growth factors in each latent class
beta_x1 <- 0.5 # coefficients for covariates X_1
means <- list(c(3,0.2), c(-3, -0.2)) # means of the growth factors in two latent classes


data_generation <- function(Nobs, 
                            Nclass, 
                            Nsubs,
                            Ninds,
                            Npops)
{
  library(dplyr)
  
  # Generate the latent class variable with the covariates
  dat <- data.frame()
  dat <- data.frame( 
    x.1 = rnorm(Npops,
                mean = 0,
                sd =1)) # add covariates
  
  for (i in 1:Npops){
    dat <- dat %>% 
      mutate(z = rbinom(Npops,
                        size =1,
                        prob = 1/(1 + exp(-(0 + beta_x1*x.1))))) 
  } # generate latent class variable influenced by covariate X.1
  
  dat_split <- split(dat, dat$z) #split two latent classes
  
  for (c in 1:Nclass){
    
    dat <- list()
    # Generated growth factors
    D <- sqrt(diag(c(1, 0.25)))                         # define the standard deviation matrix (sqrt of variances)
    R <- matrix(c(1, -0.3, -0.3, 1), nrow=2, byrow=TRUE)  # growth factors with a negative correlation of 0.3
    
    dat[[c]] <- data.frame( 
      eta = MASS::mvrnorm(nrow(dat_split[[c]]),
                          mu = means[[c]],               # mean of the growth factors
                          Sigma = D %*% R %*% D,        # covariance matrix of the factors
                          empirical = TRUE))
    
    dat_split[[c]] <- cbind(dat[[c]],dat_split[[c]])
    
    # Generated repeated measures
    for (i in (1:Ninds)){
      for (t in (1:Nobs)){
        dat_split[[c]] <- dat_split[[c]] %>% 
          mutate(y = 0 + eta.1 + t*eta.2) %>%     # generate the indicators of each time point.
          #      mutate(y = y + rnorm(Nsubs, mean=0, sd= 1/8*sd(y))) %>%  
          plyr::rename(c("y" = paste0("y", i, ".", t)))
      }
    }
  }
  
  dat_split <- bind_rows(dat_split)
  dat_split <- dat_split %>% mutate(id = 1:nrow(dat_split))   # add the id varible  
}

dat <- data_generation(Nobs, Nclass,Nsubs,Ninds, Npops)

# Save the datasets
save(dat, file = "Study1_datageneration.rda")
write.table(dat,'C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_1\\test1.csv',
            row.names=FALSE,col.names=TRUE,sep=",",quote = FALSE)

