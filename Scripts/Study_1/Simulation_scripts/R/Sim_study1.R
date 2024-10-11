# Load required package
library(dplyr)
library(brew)

# Specify the LatentGOLD syntex folder address
setwd( "~/Documents/GitHub/Two-step-GMM-code/Project_1/Scripts/Study_1/Simulation_scripts/LatentGOLD" )

# Define parameters and manipulated factors
# Regression coefficients
beta_0 <- list(c(0.75),c(1.61))
beta_x1 <- c(-0.50)
beta_x2 <- c(0.75)
beta_x1x2 <- list(c(-0.50, 0.75)) # latent class variables

# Latent class variables
alpha_0 <- list(c(3,-3), c(1.5,-1.5))
alpha_1 <- list(c(-0.2, 0.2))
gamma_0 <- list(c(0.5, -0.5))
gamma_1 <- list(c(0.25,-0.25)) # growth factors

# Variance and covariance
theta3 <- list(rep(1, 9)) # three time points
theta6 <- list(rep(1, 18)) # six time points

Sigma_1 <- list(c(1, -0.15, 1)) # for class 1
Sigma_2 <- list(c(1, -0.15, 1)) # for class 2

# Sample size
samplesize <- matrix(c('casesum500', 'casesum1000'), nrow = 2, dimnames = list(c("500", "1000")))

# Create population model parameters grid
create_parasample <- function(..., rownames_list = c("EH", "UH", "EL", "UL")) { 
  # Label the conditions: equal/unequal class proportions and high/low (relative) entropy 
  parasample <- expand.grid(...)
  rownames(parasample) <- rownames_list
  return(parasample)
}
# "EH": equal proportion  and high entropy
# "UH": unequal proportion and high entropy
# "EL": equal proportion and low(relative) entropy
# "EH": equal proportion and low(relative) entropy

# Assign parameters of population model I
parasample_pop <- create_parasample(beta_0=beta_0, beta_x1=beta_x1, 
                                    alpha_0=alpha_0, gamma_0=gamma_0, 
                                    theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)

# Assign starting values for nine estimators (prevent class switching problem)
parasample_onestep <- create_parasample(beta_0=beta_0, beta_x1=beta_x1, 
                                        alpha_0=alpha_0, alpha_1=alpha_1, 
                                        theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2) # one step models

parasample_stepone <- create_parasample(beta_0=beta_0, 
                                        alpha_0=alpha_0, alpha_1=alpha_1, 
                                        theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2) # step one models of step-wise estimators

parasample <- list(parasample_pop=parasample_pop, 
                   parasample_onestep=parasample_onestep, 
                   parasample_stepone=parasample_stepone)

# Specify the ".brew" file names of each models in different time point conditions
models_3T <- list(
  population_mode = "DatGen_3T",
  onestep = "onestep_3T",
  stepone = "step1_3T",
  steptwo = "step2_3T",
  stepthree = "step3_3T"
) # three time points

models_3T <- list(
  population_mode = "DatGen_6T",
  onestep = "onestep_6T",
  stepone = "step1_6T",
  steptwo = "step2_6T",
  stepthree = "step3_6T"
) # six time points

models <- models_3T

# Simulated data generation
library(brew) 
source("~/Documents/GitHub/Two-step-GMM-code/Project_1/Scripts/Study_1/Simulation_scripts/R/run.template.R")

for(i in 1:nsim)
{
  for(j in 1:length(samplesize))
  {
    
    for (n in 1:nrow(parasample)) {
      
      # Generating simulated data from population model I
      envir <- new.env()
      assign("samplesize", samplesize[j], envir = envir)  # label the sample size condition
      assign("sampleparameters", unlist(parasample$parasample_pop[n,]), envir = envir)  # label the entropy and proportion condition
      brewfile <- paste(population_model, ".brew", sep="")
      run.template(brewfile, envir= envir, temp.filename.base=population_model)
      lgsfile <- paste(population_model, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe", lgsfile, "/b"))
      
      
      # One-step full information maximum likelihood estimator 
      envir <- new.env()
      assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
      assign("condition", rownames(parasample)[n], envir = envir) # label the entropy and proportion condition
      assign("sampleparameters", unlist(parasample$parasample_onestep[n,]), envir = envir)  # fix the starting value at the population parameters
      brewfile1 <- paste(onestep, ".brew", sep="")
      run.template(brewfile1, envir=envir, temp.filename.base=onestep) 
      lgsfile1 <- paste(onestep, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                  lgsfile1, "/b"))	 #run LG form R 
      
      # Two-step and three-step estimators
      ## step 1 model 
      envir <- new.env()
      assign("size", rownames(samplesize)[j], envir = envir) 
      assign("sampleparameters", unlist(parasample_s1[n,]), envir = envir)  
      assign("condition", rownames(parasample_s1)[n], envir = envir) 
      brewfile2 <- paste(stepone, ".brew", sep="")
      run.template(brewfile2, envir=envir, temp.filename.base=stepone) .
      lgsfile2 <- paste(stepone, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                  lgsfile2, "/b"))
      
      # step 2 model
      DAT <-paste("step2pars",".dat", sep="")
      bpars <- as.matrix(read.table(DAT, sep=" ", dec=","))
      envir <- new.env()
      assign("condition", rownames(parasample_s1)[n], envir = envir) 
      assign("size", rownames(samplesize)[j], envir = envir) 
      assign("bpars", bpars[-1,], envir=envir)
      brewfile4 <- paste(steptwo, ".brew", sep="")
      run.template(brewfile4, envir=envir, temp.filename.base=steptwo)
      lgsfile4 <- paste(steptwo, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                  lgsfile4, "/b"))
      
      # step 3 model
      brewfile4 <- paste(stepthree, ".brew", sep="")
      run.template(brewfile4, envir=envir, temp.filename.base=stepthree).
      lgsfile4 <- paste(stepthree, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                  lgsfile4, "/b"))	
    }
  }
} 

