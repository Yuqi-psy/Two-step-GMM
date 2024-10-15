# Load required package
library(dplyr)
library(brew)

# Specify the LatentGOLD syntex folder address
setwd( )

# Define parameters and manipulated factors
  # Regression coefficients
  beta_0 <- list(-1.5, -0.43)
  beta_x1 <- c(-0.50)
  beta_x2 <- c(0.75)
  beta_x1x2 <- list(c(-0.50, 0.75)) # latent class variables
  
  # Latent class variables
  alpha_0 <- list(c(1.5, -1.5), c(0.5, -0.5))
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

# Assign parameters of population model II
parasample_pop <- create_parasample(beta_0=beta_0, beta_x1x2=beta_x1x2, 
                                     alpha_0=alpha_0, gamma_0=gamma_0, alpha_1=alpha_1, 
                                     theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)

# Assign starting values for nine estimators (prevent class switching problem)
parasample_onestepA <- create_parasample(beta_0=beta_0, beta_x1x2=beta_x1x2, 
                                         alpha_0=alpha_0, alpha_1=alpha_1, 
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)
parasample_onestepB <- create_parasample(beta_0=beta_0, beta_x1x2=beta_x1x2, 
                                         alpha_0=alpha_0, gamma_0=gamma_0, alpha_1=alpha_1, 
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)
parasample_onestepC <- create_parasample(beta_0=beta_0, beta_x1x2=beta_x1x2, 
                                         alpha_0=alpha_0, gamma_0=gamma_0, alpha_1=alpha_1, gamma_1=gamma_1,
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2) # one step models

parasample_steponeA <- create_parasample(beta_0=beta_0, 
                                         alpha_0=alpha_0, alpha_1=alpha_1, 
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)
parasample_steponeB <- create_parasample(beta_0=beta_0, beta_x2=beta_x2, 
                                         alpha_0=alpha_0, gamma_0=gamma_0, alpha_1=alpha_1, 
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2)
parasample_steponeC <- create_parasample(beta_0=beta_0, beta_x2=beta_x2, 
                                         alpha_0=alpha_0, gamma_0=gamma_0, alpha_1=alpha_1, gamma_1=gamma_1,
                                         theta=theta3, Sigma_1=Sigma_1, Sigma_2=Sigma_2) # step one models of step-wise estimators

parasample <- list(parasample_pop=parasample_pop, 
                   parasample_onestepA=parasample_onestepA, parasample_onestepB=parasample_onestepB, parasample_onestepC=parasample_onestepC, 
                   parasample_steponeA=parasample_steponeA, parasample_steponeB=parasample_steponeB, parasample_steponeC=parasample_steponeC)

# Specify the ".brew" file names of each models in different time point conditions
three_timepts <- list(
  pop_model = "DatGen_3T",
  onestep_a = "onestepA_3T", 
  onestep_b = "onestepB_3T", 
  onestep_c = "onestepC_3T", 
  steponeA = "step1A_3T", 
  steponeB = "step1B_3T", 
  steponeC = "step1C_3T", 
  steptwoA = "step2A_3T", 
  steptwoB = "step2B_3T", 
  steptwoC = "step2C_3T", 
  stepthreeA = "step3A_3T",
  stepthree = "step3_3T"
) # three time points

six_timepts <- list(
  pop_model = "DatGen_6T",
  onestep_a = "onestepA_6T", 
  onestep_b = "onestepB_6T", 
  onestep_c = "onestepC_6T", 
  steponeA = "step1A_6T", 
  steponeB = "step1B_6T", 
  steponeC = "step1C_6T", 
  steptwoA = "step2A_6T", 
  steptwoB = "step2B_6T", 
  steptwoC = "step2C_6T", 
  stepthreeA = "step3A_6T",
  stepthree = "step3_6T"
) # six time points


# Simulation
SimFnc <- function(timepoints, # conditions of time points
                   parasample, # starting values of models
                   nsim # number of simulation
) {
  library(brew)
  source(
    "/Users/liuyuqi/Documents/GitHub/Two-step-GMM-code/Project_1/Scripts/Study_2/Simulation_scripts/R/run.template.R"
  )
  # Label specification types
  modeltype <- c("A", "B", "C")
  # Specify sample size conditions
  samplesize <- matrix(c('casesum500', 'casesum1000'), nrow = 2, dimnames = list(c("500", "1000")))
  
  for (i in 1:nsim)
  {
    for (j in 1:length(samplesize))
    {
      for (n in 1:nrow(parasample$parasample_pop)) {
        
        # Generating simulated data from population model III
        envir <- new.env()
        assign("samplesize", samplesize[j], envir = envir)  # label the sample size condition
        assign("sampleparameters", unlist(parasample$parasample_pop[n, ]), envir = envir)  # label the entropy and proportion conditions
        brewfile <- paste( timepoints$pop_model, ".brew", sep = "")
        run.template(brewfile,
                     envir = envir,
                     temp.filename.base =  timepoints$pop_model)
        lgsfile <- paste( timepoints$pop_model, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile,
            "/b"
          )
        )
        
        
        # One-step full information maximum likelihood estimator (fixed starting values)
        envir <- new.env()
        assign("size", rownames(samplesize)[j], envir = envir) # label the sample size condition
        assign("condition", rownames(parasample$parasample_pop)[n], envir = envir) # label the entropy and proportion conditions
        assign("sampleparameters",
               unlist(parasample$parasample_onestepA[n, ]),
               envir = envir)  # 2) assign the population parameters
        brewfile1 <- paste(onestep_a, ".brew", sep = "")
        run.template(brewfile1,
                     envir = envir,
                     temp.filename.base = onestep_a) #create lgs files.
        lgsfile1 <- paste( timepoints$onestep_a, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile1,
            "/b"
          )
        )	 # specification A
        
        envir <- new.env()
        assign("size", rownames(samplesize)[j], envir = envir) 
        assign("condition", rownames(parasample$parasample_onestepB)[n], envir = envir) 
        assign("sampleparameters",
               unlist(parasample$parasample_onestepB[n, ]),
               envir = envir)  
        brewfile2 <- paste(timepoints$onestep_b, ".brew", sep = "")
        run.template(brewfile2,
                     envir = envir,
                     temp.filename.base = timepoints$onestep_b) #create lgs files.
        lgsfile2 <- paste( timepoints$onestep_b, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile2,
            "/b"
          )
        )	 # specification B
        
        
        envir <- new.env()
        assign("size", rownames(samplesize)[j], envir = envir)
        assign("condition", rownames(parasample$parasample_onestepC)[n], envir = envir) 
        assign("sampleparameters",
               unlist(parasample$parasample_onestepC[n, ]),
               envir = envir)  
        brewfile3 <- paste( timepoints$onestep_c, ".brew", sep = "")
        run.template(brewfile3,
                     envir = envir,
                     temp.filename.base = timepoints$onestep_c) #create lgs files.
        lgsfile3 <- paste( timepoints$onestep_c, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile3,
            "/b"
          )
        )	 # specification C
        
        # Two-step and Three-step estimators (fixed starting values)
        
        ## Specification A
        # step 1 model
        envir <- new.env()
        assign("modeltype", modeltype[1], envir = envir)
        assign("size", rownames(samplesize)[j], envir = envir) 
        assign("sampleparameters",
               unlist(parasample$parasample_steponeA[n, ]),
               envir = envir)  # 2) fix the starting values
        assign("condition", rownames(parasample$parasample_steponeA)[n], envir = envir) 
        brewfile4 <- paste( timepoints$steponeA, ".brew", sep = "")
        run.template(brewfile4,
                     envir = envir,
                     temp.filename.base = timepoints$steponeA)
        lgsfile4 <- paste( timepoints$steponeA, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile4,
            "/b"
          )
        )
        
        # step 2 model 
        DAT <- paste("step2pars_A", ".dat", sep = "")
        bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
        envir <- new.env()
        assign("condition", rownames(parasample$parasample_steponeA)[n], envir = envir)
        assign("size", rownames(samplesize)[j], envir = envir) 
        
        assign("step1a_pars", bpars[-1, ], envir = envir)
        brewfile5 <- paste( timepoints$steptwoA, ".brew", sep = "")
        run.template(brewfile5,
                     envir = envir,
                     temp.filename.base = timepoints$steptwoA) 
        lgsfile5 <- paste( timepoints$steptwoA , ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile5,
            "/b"
          )
        )	
        
        
        ## Specification B
        # step 1 model
        envir <- new.env()
        assign("modeltype", modeltype[2], envir = envir)
        assign("size", rownames(samplesize)[j], envir = envir) 
        assign("sampleparameters",
               unlist(parasample$parasample_steponeB[n, ]),
               envir = envir)  
        assign("condition", rownames(parasample$parasample_onestepB)[n], envir = envir) 
        brewfile6 <- paste( timepoints$steponeB, ".brew", sep = "")
        run.template(brewfile6,
                     envir = envir,
                     temp.filename.base = timepoints$steponeB)
        lgsfile6 <- paste( timepoints$steponeB, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile6,
            "/b"
          )
        )	 
        
        # step 2 model 
        DAT <- paste("step2pars_B", ".dat", sep = "")
        bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
        envir <- new.env()
        assign("condition", rownames(parasample$parasample_steponeB)[n], envir = envir) 
        assign("size", rownames(samplesize)[j], envir = envir) 
        assign("step1b_pars", bpars[-c(1, 2), ], envir = envir)
        
        brewfile7 <- paste( timepoints$steptwoB, ".brew", sep = "")
        run.template(brewfile7,
                     envir = envir,
                     temp.filename.base = timepoints$steptwoB) 
        lgsfile7 <- paste( timepoints$steptwoB , ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile7,
            "/b"
          )
        )	
        
        ## Specification C
        # step 1 model
        envir <- new.env()
        assign("modeltype", modeltype[3], envir = envir)
        assign("size", rownames(samplesize)[j], envir = envir)
        assign("sampleparameters",
               unlist(parasample$parasample_steponeC[n, ]),
               envir = envir)  
        assign("condition", rownames(parasample$parasample_steponeC)[n], envir = envir) 
        brewfile8 <- paste( timepoints$steponeC, ".brew", sep = "")
        run.template(brewfile8,
                     envir = envir,
                     temp.filename.base = timepoints$steponeC)
        lgsfile8 <- paste( timepoints$steponeC, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile8,
            "/b"
          )
        )=
          
          # step 2 model
          DAT <- paste("step2pars_C", ".dat", sep = "")
        bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
        envir <- new.env()
        assign("condition", rownames(parasample$parasample_steponeC)[n], envir = envir) 
        assign("size", rownames(samplesize)[j], envir = envir) 
        assign("step1c_pars", bpars[-c(1, 2), ], envir = envir)
        
        ## step 2 B model
        brewfile9 <- paste( timepoints$steptwoC, ".brew", sep = "")
        run.template(brewfile9,
                     envir = envir,
                     temp.filename.base = timepoints$steptwoC)
        lgsfile9 <- paste( timepoints$steptwoC , ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile9,
            "/b"
          )
        )	
        
        # step 3 model
        envir <- new.env()
        assign("modeltype", modeltype[1], envir = envir)
        assign("stepone", modeltype[1], envir = envir)
        assign("condition",
               rownames(parasample$parasample_steponeA)[n],
               envir = envir)
        assign("size", rownames(samplesize)[j], envir = envir)
        brewfile10 <- paste( parasample$timepoints$stepthreeA, ".brew", sep = "")
        run.template(brewfile10,
                     envir = envir,
                     temp.filename.base = timepoints$stepthreeA) 
        lgsfile10 <- paste( timepoints$stepthreeA, ".lgs", sep = "")
        shell(
          paste(
            "C:~\\lg60.exe",
            lgsfile10,
            "/b"
          )
        ) # Specification A
        
        for (s1_model in modeltype[-1]) {
          envir <- new.env()
          assign("modeltype", s1_model, envir = envir)
          assign("stepone", s1_model, envir = envir)
          assign("condition",
                 rownames(parasample$parasample_steponeA)[n],
                 envir = envir)
          assign("size", rownames(samplesize)[j], envir = envir)
          brewfile <- paste( timepoints$stepthree, ".brew", sep = "")
          run.template(brewfile,
                       envir = envir,
                       temp.filename.base = timepoints$stepthree) 
          lgsfile <- paste( timepoints$stepthree, ".lgs", sep = "")
          shell(
            paste(
              "C:~\\lg60.exe",
              lgsfile,
              "/b"
            )
          )	 
        } # Specification B and C
      }
    }
  }
}

# Running simulation for study 2 at 3 time point condition
RunSim <- SimFnc(timepoints = three_timepts, parasample = parasample, nsim = 100)



