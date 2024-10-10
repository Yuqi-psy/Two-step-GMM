library(dplyr)

########## Study 2 population model parameters ##########

Int_lc <- list(c(-1.5),c(-0.43)) # the intercept of the latent class variable for equal(-1.9) and unequal(-0.80) proportion
Slo_lc_cov1 <- c(-0.50)
Slo_lc_cov12 <- list(c(-0.50, 0.75)) # the slope of latent class variable, i.e., covariates effects of X1(fixed)
Slo_lc_cov2 <- c(0.75)

Int_i0 <- list(c(1.5,-1.5), c(0.5,-0.5)) # the intercept of latent intercept factor for high(3,-3) and low(1,-1) entropy
Int_i1 <- list(c(0.5,-0.5))
Slo_s0 <- list(c(-0.2,0.2)) # the slope of latent slope factor
Slo_s1 <- list(c(0.25,-0.25)) # the coefficients for cov2 on the latent slope

Var_items3 <- list(c(rep(1,9))) # variance of items for 3 time points
Var_items6 <- list(c(rep(1,18))) # variance of items for 6 time points

Sigma_1 <- list(c(1, -0.15, 1)) # variance-covariance of growth factors for latent class 1
Sigma_2 <- list(c(1, -0.15, 1)) # variance-cpvariance of growth factors for latent class 2

conditions <- c("EH", # "EH": equal class proportions  and high entropy
                "UH", # "UH": unequal class proportions and high entropy
                "EL", # "EL": equal class proportions and low entropy
                "UL") # "EH": equal class proportions and low entropy




# Create population model parameters grid
parasample_pop <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc_cov12, # regression coefficient for the latent class variable
                              Int_i0=Int_i0, Int_i1=Int_i1, # regression coefficient for the latent intercept variable
                              Slo_s0=Slo_s0, # regression coefficient for the latent slope variable
                               Var_items=Var_items6,
                              Sigma_1=Sigma_1, Sigma_2=Sigma_2) %>% rownames(conditions)
rownames(parasample_pop) <- c("EH", "UH", "EL", "UL") # "EH": equal proportion  and high entropy
                                                      # "UH": unequal proportion and high entropy
                                                      # "EL": equal proportion and low entropy
                                                      # "EH": equal proportion and low entropy
# Create population model parameters grid (step one model)

parasample_onestepA <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc_cov12,
                                   Int_li=Int_i0, 
                                   Slo_s0=Slo_s0,
                                   Var_items=Var_items6,
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2
                                   # !!Remember to specify the time points
)
rownames(parasample_onestepA) <- c("EH", "UH", "EL", "UL") 

parasample_onestepB <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc_cov12,
                                   Int_li=Int_i0, Int_i1=Int_i1, 
                                   Slo_s0=Slo_s0, 
                                   Var_items=Var_items6,
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2
                                   # !!Remember to specify the time points
)
rownames(parasample_onestepB) <- c("EH", "UH", "EL", "UL") 

parasample_onestepC <- expand.grid(Int_lc=Int_lc, Slo_lc=Slo_lc_cov12,
                                   Int_li=Int_i0, Int_i1=Int_i1, 
                                   Slo_s0=Slo_s0, Slo_s1=Slo_s1,
                                   Var_items=Var_items6,
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2
                                   # !!Remember to specify the time points
                                  )
rownames(parasample_onestepC) <- c("EH", "UH", "EL", "UL") 

parasample_steponeA <- expand.grid(Int_lc=Int_lc,
                                   Int_i0=Int_i0,
                                   Slo_s0=Slo_s0, 
                                   Var_items=Var_items6,
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2# !!Remember to specify the time points
                                   )
rownames(parasample_steponeA) <- c("EH", "UH", "EL", "UL")

parasample_steponeB <- expand.grid(Int_lc=Int_lc, Slo_lc_cov1=Slo_lc_cov2,
                                   Int_i0=Int_i0, Int_i1=Int_i1,
                                   Slo_s0=Slo_s0,
                                   Var_items=Var_items6,  
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2# !!Remember to specify the time points
                                   )
rownames(parasample_steponeB) <- c("EH", "UH", "EL", "UL") 

parasample_steponec <- expand.grid(Int_lc=Int_lc, Slo_lc_cov1=Slo_lc_cov2,
                                   Int_i0=Int_i0, Int_i1=Int_i1,
                                   Slo_s0=Slo_s0, Slo_s1=Slo_s1,
                                   Var_items=Var_items6, 
                                   Sigma_1=Sigma_1, Sigma_2=Sigma_2# !!Remember to specify the time points
                                   )
rownames(parasample_steponec) <- c("EH", "UH", "EL", "UL") 

parasample <- list(parasample_pop=parasample_pop, parasample_onestepA=parasample_onestepA, parasample_onestepB = parasample_onestepB, parasample_onestepC=parasample_onestepC,
                   parasample_steponeA=parasample_steponeA, parasample_steponeB=parasample_steponeB, parasample_steponec=parasample_steponec)

samplesize <- matrix(c('casesum500', 'casesum1000')) # the total sample size for population model
rownames(samplesize) <- c("500","1000")

########### GMM_simulation ###########

# set the model working directory
setwd("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_2\\6T")


# specify models
population_model <- c("DatGen_6T") # population model

onestep_a <- c("onestepA_6T")
onestep_b <- c("onestepB_6T")
onestep_c <- c("onestepC_6T") # one step model

steponeA <- c("step1A_6T")
steponeB <- c("step1B_6T") 
steponeC <- c("step1C_6T") # step one model

steptwoA1 <- c("step2A1_6T")
steptwoA2 <- c("step2A2_6T") 
steptwoB <- c("step2B_6T") # step two model

stepthree <- c("step3_6T") # step three model

modeltype <- c("A", "B", "C") # step one model type

# define number of simulaions
nsim <- 99

######

for(i in 1:nsim)
     {
  st <- Sys.time()
   for(j in 1:length(samplesize))
      { # run three variance condition in one sim
     library(brew) 
     source("C:\\Users\\liuy16\\Desktop\\PhD files\\Project1\\Simulated_data\\Study_1\\3_time_points\\run.template.R")
     
  for (n in 1:nrow(parasample$parasample_pop)) {

    
  # Data generation
    envir <- new.env()
    assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("sampleparameters", unlist(parasample$parasample_pop[n, ]), envir = envir)  # 2) assign the population parameters
    brewfile <- paste(population_model, ".brew", sep="")
    run.template(brewfile, envir= envir, temp.filename.base=population_model)
    lgsfile <- paste(population_model, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe", lgsfile, "/b"))
    
    
  # One-step full information ML model (fixed starting values)
    # model A
    envir <- new.env()
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("condition", rownames(parasample$parasample_pop)[n], envir = envir) # label the entropy and proportion condition
    assign("sampleparameters", unlist(parasample$parasample_onestepA[n,]), envir = envir)  # 2) assign the population parameters
    brewfile1 <- paste(onestep_a, ".brew", sep="")
    run.template(brewfile1, envir=envir, temp.filename.base=onestep_a) #create lgs files.
    lgsfile1 <- paste(onestep_a, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile1, "/b"))	 #run LG form R 
    # !!!!For assign sample parameter problem, create list of sample parameter grids for different population model and step one model
    # model B
    envir <- new.env()
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("condition", rownames(parasample$parasample_pop)[n], envir = envir) # label the entropy and proportion condition
    assign("sampleparameters", unlist(parasample$parasample_onestepB[n,]), envir = envir)  # 2) assign the population parameters
    brewfile1 <- paste(onestep_b, ".brew", sep="")
    run.template(brewfile1, envir=envir, temp.filename.base=onestep_b) #create lgs files.
    lgsfile1 <- paste(onestep_b, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile1, "/b"))	 #run LG form R 
    
    # model C
    envir <- new.env()
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("condition", rownames(parasample$parasample_pop)[n], envir = envir) # label the entropy and proportion condition
    assign("sampleparameters", unlist(parasample$parasample_onestepC[n,]), envir = envir)  # 2) assign the population parameters
    brewfile1 <- paste(onestep_c, ".brew", sep="")
    run.template(brewfile1, envir=envir, temp.filename.base=onestep_c) #create lgs files.
    lgsfile1 <- paste(onestep_c, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile1, "/b"))	 #run LG form R

  # Two-step approach
    ##### METHOD 1
    # step 1 model A
    envir <- new.env()
    #assign("modeltype_s2", s2_A1A2, envir = envir)
    assign("modeltype", modeltype[1], envir = envir)
    #assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("sampleparameters", unlist(parasample$parasample_steponeA[n,]), envir = envir)  # 2) fix the starting values
    assign("condition", rownames(parasample$parasample_steponeA)[n], envir = envir) # label the entropy and proportion condition
    brewfile2 <- paste(steponeA, ".brew", sep="")
    run.template(brewfile2, envir=envir, temp.filename.base=steponeA) #create lgs files.
    lgsfile2 <- paste(steponeA, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile2, "/b"))	 #run LG form R 
    
    # step 2 model (Step one model: step one A)
    DAT <-paste("step2pars_A",".dat", sep="")
    bpars <- as.matrix(read.table(DAT, sep=" ", dec=","))
    envir <- new.env()
    assign("condition", rownames(parasample$parasample_steponeA)[n], envir = envir) # label the entropy and proportion condition
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    
    
    ## step 2 A1 model(ignore the DF)
    assign("step1a_pars", bpars[-1, ], envir=envir)
    brewfile4 <- paste(steptwoA1, ".brew", sep="")
    run.template(brewfile4, envir=envir, temp.filename.base=steptwoA1) #create lgs files.
    lgsfile4 <- paste(steptwoA1 , ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile4, "/b"))	 #run LG form R 
    
    
    ##### METHOD 2
    # step 1 model B
    envir <- new.env()
    #assign("modeltype_s2", s2_A1, envir = envir)
    assign("modeltype", modeltype[2], envir = envir)
    #assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("sampleparameters", unlist(parasample$parasample_steponeB[n,]), envir = envir)  # 2) fix the starting values
    assign("condition", rownames(parasample$parasample_onestepB)[n], envir = envir) # label the entropy and proportion condition
    brewfile2 <- paste(steponeB, ".brew", sep="")
    run.template(brewfile2, envir=envir, temp.filename.base=steponeB) #create lgs files.
    lgsfile2 <- paste(steponeB, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile2, "/b"))	 #run LG form R 
    
    # step 2 model (Step one model: step one B)
    DAT <-paste("step2pars_B",".dat", sep="")
    bpars <- as.matrix(read.table(DAT, sep=" ", dec=","))
    envir <- new.env()
    assign("condition", rownames(parasample$parasample_steponeB)[n], envir = envir) # label the entropy and proportion condition
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("step1b_pars", bpars[-c(1,2),], envir=envir)
    
    ## step 2 A1 model(ignore the DF)
    brewfile4 <- paste(steptwoA2, ".brew", sep="")
    run.template(brewfile4, envir=envir, temp.filename.base=steptwoA2) #create lgs files.
    lgsfile4 <- paste(steptwoA2 , ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile4, "/b"))	 #run LG form R 
  
    ##### METHOD 3
    # step 1 model C
    envir <- new.env()    
    #assign("modeltype_s2", s2_B, envir = envir)
    assign("modeltype", modeltype[3], envir = envir)
    #assign("samplesize", samplesize[j], envir = envir)  # 1) assign the condition of sample sizes
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("sampleparameters", unlist(parasample$parasample_steponec[n,]), envir = envir)  # 2) fix the starting values
    assign("condition", rownames(parasample$parasample_steponec)[n], envir = envir) # label the entropy and proportion condition
    brewfile2 <- paste(steponeC, ".brew", sep="")
    run.template(brewfile2, envir=envir, temp.filename.base=steponeC) #create lgs files.
    lgsfile2 <- paste(steponeC, ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile2, "/b"))	 #run LG form R
    
    # step 2 model (Step one model: step one C)
    DAT <-paste("step2pars_C",".dat", sep="")
    bpars <- as.matrix(read.table(DAT, sep=" ", dec=","))
    envir <- new.env()
    assign("condition", rownames(parasample$parasample_steponec)[n], envir = envir) # label the entropy and proportion condition
    assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
    assign("step1c_pars", bpars[-c(1,2),], envir=envir)
    
    ## step 2 B model
    brewfile4 <- paste(steptwoB, ".brew", sep="")
    run.template(brewfile4, envir=envir, temp.filename.base=steptwoB) #create lgs files.
    lgsfile4 <- paste(steptwoB , ".lgs", sep="")
    shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                lgsfile4, "/b"))	 #run LG form R 
 
    # step 3 model
    for (s1_model in modeltype) {
      print(s1_model)
      envir <- new.env()    
      assign("modeltype", s1_model, envir = envir)
      assign("stepone", s1_model, envir = envir)
      assign("condition", rownames(parasample$parasample_steponeA)[n], envir = envir)
      assign("size", rownames(samplesize)[j], envir = envir)
      brewfile5 <- paste(stepthree, ".brew", sep="")
      run.template(brewfile5, envir=envir, temp.filename.base=stepthree) #create lgs files.
      lgsfile4 <- paste(stepthree, ".lgs", sep="")
      shell(paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe",
                  lgsfile4, "/b"))	 #run LG form R 
     }
    }
   }
  et <- Sys.time()
  } 

st-et
