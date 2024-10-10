



source("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\R_latentgold\\run.template.R") # using this source file you create the .lgs files that will be run in LG. using this option is easy to create various conditions usinga single file
source("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\latentgold\\data_generation.R")
library(brew) # need to generate .lgs files 

########## Study 1 population model parameters ##########
samplesize <- matrix(c('casesum500', 'casesum1000')) # the total sample size for population model
rownames(samplesize) <- c("500","1000")

Int_lc <- list(c(0.75),c(1.61)) # the intercept of the latent class variable for equal(0.5) and unequal(0.32) proportion
Slo_lc <- c(-0.5) # the slope of latent class variable, i.e., covariates effects of X1(fixed)

Int_li <- list(c(3,-3), c(1.5,-1.5)) # the intercept of latent intercept factor for high(3,-3) and low(1,-1) entropy
Slo_ls <- list(c(-0.2,0.2)) # the slope of latent slope factor

Var_items3 <- list(c(rep(1,9))) # variance of items for 3 time points
Var_items6 <- list(c(rep(1,18))) # variance of items for 6 time points

Sigma_1 <- list(c(1, -0.15, 1)) # variance-covariance of growth factors for latent class 1
Sigma_2 <- list(c(1, -0.15, 1)) # variance-cpvariance of growth factors for latent class 2

# Create population model parameters grid 
parasample <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc, Int_li=Int_li, 
                          Slo_ls=Slo_ls, Var_items=Var_items6, # !!Remember to specify the time points
                          Sigma_1=Sigma_1, Sigma_2=Sigma_2)
rownames(parasample) <- c("EH", "UH", "EL", "UL") # "EH": equal proportion  and high entropy
                                                  # "UH": unequal proportion and high entropy
                                                  # "EL": equal proportion and low entropy
                                                  # "EH": equal proportion and low entropy
# Create population model parameters grid (step one model)
parasample_s1 <- expand.grid(Int_lc=Int_lc, Int_li=Int_li, 
                          Slo_ls=Slo_ls, Var_items=Var_items6, # !!Remember to specify the time points
                          Sigma_1=Sigma_1, Sigma_2=Sigma_2)
rownames(parasample_s1) <- c("EH", "UH", "EL", "UL") # "EH": equal proportion  and high entropy
                                                     # "UH": unequal proportion and high entropy
                                                     # "EL": equal proportion and low entropy
                                                     # "EL": equal proportion and low entropy

########################loop over time coding using brew file + sim condition of sampelsize #####################




lc <- c(0.5, -0.75) # coef of latent class
covariate_c <- c(-1,1) # coef for covariate in latent class
covatiate_s <- c(-1.5,0,1.5) # coef for covariate in latent slope

covariance <- c(-1.0,0,1.0) # covariance for growth factors
intercept <- c(-1.5,0,1.5) # coef for latent intercept
slope <- c(-0.5,-1.0,1.0) # coef for latent slope

parameters <-c(lc, covariate_c, variance, covariance, intercept, slope, covariate_s)  ## create the matrix of parameters

# define population model parameters
var_low <- c(0,0,0,0,0,0)
var_middle <- c(1,1,0,1,1,0)
var_high <- c(1,1,2,1,1,2) # question need to be discuss, the variance of which class should be fix to 0? in terms of what
var_condition <- cbind(var_low,var_middle,var_high)

# define the population model
pop_models <- c("GMM_gen_2cov", "GMM_gen_covs", "GMM_gen_covc")

# define the statisitcal model
onestep_models <- c("GMM_2cov_onestep", "GMM_covs_onestep", "GMM_covc_onestep")
stepone_models <- c("GMM_2cov_step1", "GMM_covs_step1", "GMM_covc_step1")
steptwo_models <- c("GMM_2cov_step2", "GMM_covs_step2", "GMM_covc_step2")

# define the fixed parameters
cov_lc_label <- c("c d", "c d e", "c d e f", "c d e f g" )
cov_s_label <- c("a", "a c d", "c d", "c d e", "c d e f","c f")
cov2_label <- c("c d", "c d e", "c d e f", "c d e f g", "c f")
cov_lc_par <- list(c(var_i, var_s), c(var_i, var_s, covar), c(var_i, var_s,covar,coef_int),c(var_i, var_s,covar,coef_int,coef_slo))
cov_s_par <- list(c(lc), c(lc,var_i, var_s), c(var_i, var_s), c(var_i, var_s,covar), c(var_i, var_s,covar,coef_int), c(var_i, coef_int))
cov2_par <- list(c(var_i, var_s), c(var_i, var_s, covar), c(var_i, var_s,covar,coef_int),c(var_i, var_s,covar,coef_int,coef_slo), c(var_i, coef_int))

select_labels <- c("c f", "c f", "c d e f g")
select_models <- list(c(var_i,coef_int), c(var_i,coef_int), c(var_i, var_s, covar,coef_int,coef_slo)) # selected optimal models for 2cov, covs,covc
########### GMM_simulation ###########

# set the model working directory
setwd("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_1\\1.5_6")
# set the number of simulation and models(population model, one step, two step)

nsim <- 98


# Default
##  Create population model parameters grid 
parasample <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc, Int_li=Int_li, 
                          Slo_ls=Slo_ls, Var_items=Var_items6, # !!Remember to specify the time points
                          Sigma_1=Sigma_1, Sigma_2=Sigma_2)
rownames(parasample) <- c("EH", "UH", "EL", "UL") 
# "EH": equal proportion  and high entropy
# "UH": unequal proportion and high entropy
# "EL": equal proportion and low entropy
# "UL": unequal proportion and low entropy

# Create step 1 model parameters grid (for fixing starting value, avoiding class switching)
parasample_s1 <- expand.grid(Int_lc=Int_lc, Int_li=Int_li, 
                             Slo_ls=Slo_ls, Var_items=Var_items6, # !!Remember to specify the time points
                             Sigma_1=Sigma_1, Sigma_2=Sigma_2)
rownames(parasample_s1) <- c("EH", "UH", "EL", "UL") 



##  specify models
population_model <- c("DatGen_6T")
onestep <- c("onestep_6T")
stepone <- c("step1_6T")
steptwo <- c("step2_6T")
stepthree <- c("step3_6T")



# Simulated data generation
library(brew) 
source("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_1\\3_time_points\\run.template.R")

######
for(i in 1:nsim)
     {
 st <- Sys.time()
   for(j in 1:length(samplesize))
      { # run three variance condition in one sim
     library(brew) 
     source("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_1\\3_time_points\\run.template.R")
     
  for (n in 1:nrow(parasample)) {

    
  # Data generation
    envir <- new.env()
    assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("sampleparameters", unlist(parasample[n,]), envir = envir)  # 2) assign the population parameters
    brewfile <- paste(population_model, ".brew", sep="")
    run.template(brewfile, envir= envir, temp.filename.base=population_model)
    lgsfile <- paste(population_model, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe", lgsfile, "/b"))
    
    
  # One-step full information ML model (fixed starting values)
    envir <- new.env()
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("condition", rownames(parasample)[n], envir = envir) # label the entropy and proportion condition
    assign("sampleparameters", unlist(parasample[n,]), envir = envir)  # 2) assign the population parameters
    brewfile1 <- paste(onestep, ".brew", sep="")
    run.template(brewfile1, envir=envir, temp.filename.base=onestep) #create lgs files.
    lgsfile1 <- paste(onestep, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile1, "/b"))	 #run LG form R 
    
  # Two-step and three-step approach
    # step 1 model 
    envir <- new.env()
    #assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("sampleparameters", unlist(parasample_s1[n,]), envir = envir)  # 2) fix the starting values
    assign("condition", rownames(parasample_s1)[n], envir = envir) # label the entropy and proportion condition
    brewfile2 <- paste(stepone, ".brew", sep="")
    run.template(brewfile2, envir=envir, temp.filename.base=stepone) #create lgs files.
    lgsfile2 <- paste(stepone, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile2, "/b"))	 #run LG form R 
    
    # step 2 model
    DAT <-paste("step2pars",".dat", sep="")
    bpars <- as.matrix(read.table(DAT, sep=" ", dec=","))
    envir <- new.env()
    assign("condition", rownames(parasample_s1)[n], envir = envir) # label the entropy and proportion condition
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("bpars", bpars[-1,], envir=envir)
    brewfile4 <- paste(steptwo, ".brew", sep="")
    run.template(brewfile4, envir=envir, temp.filename.base=steptwo) #create lgs files.
    lgsfile4 <- paste(steptwo, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile4, "/b"))	 #run LG form R 
    
    # step 3 model
    brewfile4 <- paste(stepthree, ".brew", sep="")
    run.template(brewfile4, envir=envir, temp.filename.base=stepthree) #create lgs files.
    lgsfile4 <- paste(stepthree, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile4, "/b"))	 #run LG form R 
    }
   }
 et <-Sys.time()
  } 

st-et




##### draft ####
# 1. 前景提要：忘记指定参数顺序；三步法没有标准误的结果；SE对于前类别是一样的为什么？；估计的项目方差是残差吗？
# 2. 结果查看：我们需要汇报其他的参数吗？（潜类别的参数恢复等）；有关type I error 如何衡量？ERR？EDR？（packages）
# 3. 下周计划：完成study2 的参数指定（1.协变量效应的参数确定以及开始运行3个时间点的条件）











 