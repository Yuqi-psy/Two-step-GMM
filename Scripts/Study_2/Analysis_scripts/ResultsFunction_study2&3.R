##########
library(knitr)
library(kableExtra)
library(stargazer)
library(SimDesign)
library(DescTools)
setwd( "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results")

#### default settings ####
sim_conditions <- c("EH1000", "EH500", 
                    "EL1000", "EL500",
                    "UH1000", "UH500",
                    "UL1000", "UL500")
sim_models <- c("One step A", "One step B", "One step C")

#### folder path(win) ####

folder_6T2 <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\onestepA",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\onestepB",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\onestepC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step1B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step1C",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step2A1A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step2A2B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step2BC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step3_A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step3_B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\6T\\step3_C"
                        )) # study 2

folder_3T2 <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\onestepA",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\onestepB",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\onestepC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step1B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step1C",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step2A1A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step2A2B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step2BC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step3_A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step3_B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study2\\3T\\step3_C"
                        )) # study 2

folder_6T3 <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\onestepA",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\onestepB",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\onestepC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step1B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step1C",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step2A1A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step2A2B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step2BC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step3_A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step3_B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\6T\\step3_C"
                        )) # study 3

folder_3T3 <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\onestepA",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\onestepB",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\onestepC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step1B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step1C",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step2A1A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step2A2B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step2BC",
                          
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step3_A",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step3_B",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study3\\3T\\step3_C"
                        )) #study 3


#### create folder path(mac) ####
# Study 2
folder_6T <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/onestepA",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/onestepB",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/onestepC",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step1B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step1C",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step2A1A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step2A2B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step2BC",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step3_A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step3_B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/6T/step3_C"
                        ))

folder_3T <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/onestepA",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/onestepB",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/onestepC",
                        
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step1B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step1C",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step2A1A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step2A2B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step2BC",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step3_A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step3_B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study2/3T/step3_C"
                        ))

# Study 3
folder_6T <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                               "Step1B", "Step1C",
                               "Step2A1A", "Step2A2B", "Step2BC",
                               "Step3_A","Step3_B","Step3_C"),
                     Path = c(
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/onestepA",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/onestepB",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/onestepC",
                       
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step1B",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step1C",
                       
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step2A1A",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step2A2B",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step2BC",
                       
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step3_A",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step3_B",
                       "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/6T/step3_C"
                       ))

folder_3T <- data.frame(Model = c("onestepA", "onestepB", "onestepC", 
                                  "Step1B", "Step1C",
                                  "Step2A1A", "Step2A2B", "Step2BC",
                                  "Step3_A","Step3_B","Step3_C"),
                        Path = c(
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/onestepA",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/onestepB",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/onestepC",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step1B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step1C",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step2A1A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step2A2B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step2BC",
                          
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step3_A",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step3_B",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study3/3T/step3_C"
                        ))

#### read file and complie results of each model ####

# define the folder path
folder <- folder_6T


  for (i in 1:nrow(folder)) {
  print(i)
  
    folder_path <-folder[i,2]
    # 读取文件夹中所有CSV文件
    results <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE) %>% 
      lapply( read.csv, header = F, sep = ",")
    names(results) <- list.files(folder_path, pattern = "\\.csv$", full.names = F)
    
    model <- folder[i,1]
    
    ##### select data #####
    if (model== "onestepA"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      

    
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
      } 
      
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      OnestepA_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      onestepa_cvg <- CR
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      OnestepA_cr <- CR
    
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      OnestepA_mse <- MSE
    
    }
    
    if (model== "onestepB"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Entropy <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        # Entropy
        Entropy[ ,i] <- results[[i]][c(seq(3, nrow(results[[ i ]]), by = 10)), 10]
        
        # Parameters
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        # The coefficient of latent intercept 
        Beta_0i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 9]
        Beta_0i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 10]
        
        # The coefficient of latent intercept 
        Beta_1i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 11]
        Beta_1i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 12]
        
        # The coefficient for latent slope
        Beta_0s$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 13]
        Beta_0s$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 14]
        
        
        # The residual variance for iterms
        for (j in 1:c(timepoints*3)){
          theta_3[[j]][ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 14+j]
        }
        
        
        # The variance and covariance matrix for growth factors
        #### class 1
        Sigma_1$int[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 24]
        Sigma_1$int_slo[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 25]
        Sigma_1$slo[ ,i] <-  results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 26]
        
        #### class 2
        Sigma_2$int[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 27]
        Sigma_2$int_slo[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 28]
        Sigma_2$slo[ ,i] <-  results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 29]
        
        # standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X1)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
        # The coefficient of latent intercept 
        SE_1i$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 11]
        SE_1i$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 12]
        
        
        
        
        
      } 
      ## Beta_1c and Beta_2c
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      OnestepB_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      onestepb_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      OnestepB_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      OnestepB_mse <- MSE
      
      ## Beta_1i
      # Bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1i$class1, parameter = 0.5, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1i$class1, parameter = 0.5, abs = T),
        RSE_c1 = RSE(SE_1i$class1, Beta_1i$class1),
        PRB_c2 = bias(Beta_1i$class2, parameter = -0.5, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1i$class2, parameter = -0.5, abs = T),
        RSE_c2 = RSE(SE_1i$class2, Beta_1i$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      OnestepB_1i <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1i$class1[i,j] - 1.96*SE_1i$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1i$class1[i,j] + 1.96*SE_1i$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1i$class2[i,j] - 1.96*SE_1i$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1i$class2[i,j] + 1.96*SE_1i$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.5))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      onestepB_1i_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1i$class1, 0.5, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1i$class2, -0.5, MSE = T))
      colnames(MSE) <- sim_conditions 
      onestepB_1i_mse <- MSE
    }
    
    if (model== "onestepC"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      Entropy_s3 <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        # Entropy
        Entropy_s3[ ,i] <- results[[i]][c(seq(3, nrow(results[[ i ]]), by = 10)), 10]
        
        
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        # The coefficient of latent intercept 
        Beta_0i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 9]
        Beta_0i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 10]
        
        # The coefficient of latent intercept 
        Beta_1i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 11]
        Beta_1i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 12]
        
        # The coefficient for latent slope
        Beta_0s$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 13]
        Beta_0s$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 14]
        
        # The coefficient of latent intercept 
        Beta_1s$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 15]
        Beta_1s$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 16]  
        
        # The residual variance for iterms
        for (j in 1:c(timepoints*3)){
          theta_3[[j]][ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 16+j]
          
          
          # Wald statistics
          Wald_1s[ ,i] <- results[[i]][c(seq(4, nrow(results[[ i ]]), by = 10)), 15]
          
        }
        
        
        # The variance and covariance matrix for growth factors
        #### class 1
        Sigma_1$int[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 26]
        Sigma_1$int_slo[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 27]
        Sigma_1$slo[ ,i] <-  results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 28]
        
        #### class 2
        Sigma_2$int[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 29]
        Sigma_2$int_slo[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 30]
        Sigma_2$slo[ ,i] <-  results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 31]
        
        # standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X1)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
        # The coefficient of latent intercept 
        SE_1i$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 11]
        SE_1i$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 12]
        
        # The coefficient for latent slope
        SE_1s$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 15]
        SE_1s$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 16]
      } 
      ## Beta_1c and Beta_2c
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      OnestepC_bias <- onestepB_Beta_c 
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      onestepc_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      OnestepC_cr <- CR
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      OnestepC_mse <- MSE
      
      ## Beta 1s
      # type I
      rates <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
      for (j in 1:8) {
        for (i in 1:100) {
          ifelse(Wald_1s[i,j] >= 5.99, rates[i,j] <- 1, rates[i,j] <- 0)
        }
      }
      colnames(rates) <- sim_conditions
      typeI_onestepC <- colSums(rates)
      
      ## Beta 1i
      # Bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1i$class1, parameter = 0.5, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1i$class1, parameter = 0.5, abs = T),
        RSE_c1 = RSE(SE_1i$class1, Beta_1i$class1),
        PRB_c2 = bias(Beta_1i$class2, parameter = -0.5, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1i$class2, parameter = -0.5, abs = T),
        RSE_c2 = RSE(SE_1i$class2, Beta_1i$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      OnestepC_1i <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1i$class1[i,j] - 1.96*SE_1i$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1i$class1[i,j] + 1.96*SE_1i$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1i$class2[i,j] - 1.96*SE_1i$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1i$class2[i,j] + 1.96*SE_1i$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.5))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      onestepC_1i_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1i$class1, 0.5, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1i$class2, -0.5, MSE = T))
      colnames(MSE) <- sim_conditions
      onestepC_1i_mse <- MSE
      
      ## Beta 1s
      # Bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1s$class1, parameter = 0.25, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1s$class1, parameter = 0.25, abs = T),
        RSE_c1 = RSE(SE_1s$class1, Beta_1s$class1),
        PRB_c2 = bias(Beta_1s$class2, parameter = -0.25, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1s$class2, parameter = -0.25, abs = T),
        RSE_c2 = RSE(SE_1s$class2, Beta_1s$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      OnestepC_1s <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1s$class1[i,j] - 1.96*SE_1s$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1s$class1[i,j] + 1.96*SE_1s$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1s$class2[i,j] - 1.96*SE_1s$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1s$class2[i,j] + 1.96*SE_1s$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.25))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.25))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      onestepC_1s_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1s$class1, 0.25, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1s$class2, -0.25, MSE = T))
      colnames(MSE) <- sim_conditions
      onestepC_1s_mse <- MSE
    } 
    
    if (model== "Step2A1A"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
      } 
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step2A1_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      steptwoa1_cvg <-CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      SteptwoA1_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      SteptwoA1_mse <- MSE
    }
    
    if (model== "Step1B"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)) {
        ### paramters
        # The coefficient of latent intercept 
        Beta_1i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 10]
        Beta_1i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 11]
        
        ### Standard error
        SE_1i$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 10]
        SE_1i$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 11]  
        
        
      }
      ## Beta_1i
      # bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1i$class1, parameter = 0.5, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1i$class1, parameter = 0.5, abs = T),
        RSE_c1 = RSE(SE_1i$class1, Beta_1i$class1),
        PRB_c2 = bias(Beta_1i$class2, parameter = -0.5, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1i$class2, parameter = -0.5, abs = T),
        RSE_c2 = RSE(SE_1i$class2, Beta_1i$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      step1B_1i <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1i$class1[i,j] - 1.96*SE_1i$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1i$class1[i,j] + 1.96*SE_1i$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1i$class2[i,j] - 1.96*SE_1i$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1i$class2[i,j] + 1.96*SE_1i$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.5))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      steponeB_1i_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1i$class1, 0.5, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1i$class2, -0.5, MSE = T))
      colnames(MSE) <- sim_conditions
      steponeB_1i_mse <- MSE
    }
    
    if (model== "Step2A2B"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
        
      }
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step2A2B_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      step2a2b_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      step2A2B_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      step2A2B_mse <- MSE
    }
    
    if (model== "Step1C"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        # Parameter
        Beta_1i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 10]
        Beta_1i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 11]
        
        Beta_1s$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 14]
        Beta_1s$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 15]
        
        # Standard error
        SE_1i$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 10]
        SE_1i$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 11]
        
        SE_1s$class1[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 14]
        SE_1s$class2[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 15]
        ###### One-step A
        
        # Wald statistics
        Wald_1s[ ,i] <- results[[i]][c(seq(4, nrow(results[[ i ]]), by = 10)), 14]
        
      }
      ## Beta_1i
      # bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1i$class1, parameter = 0.5, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1i$class1, parameter = 0.5, abs = T),
        RSE_c1 = RSE(SE_1i$class1, Beta_1i$class1),
        PRB_c2 = bias(Beta_1i$class2, parameter = -0.5, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1i$class2, parameter = -0.5, abs = T),
        RSE_c2 = RSE(SE_1i$class2, Beta_1i$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      step1C_1i <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1i$class1[i,j] - 1.96*SE_1i$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1i$class1[i,j] + 1.96*SE_1i$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1i$class2[i,j] - 1.96*SE_1i$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1i$class2[i,j] + 1.96*SE_1i$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.5))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      steponeC_1i_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1i$class1, 0.5, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1i$class2, -0.5, MSE = T))
      colnames(MSE) <- sim_conditions
      steponeC_1i_mse <- MSE
      
      ## Beta 1s
      # type I
      rates <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
      for (j in 1:8) {
        for (i in 1:100) {
          ifelse(Wald_1s[i,j] >= 5.99, rates[i,j] <- 1, rates[i,j] <- 0)
        }
      }
      colnames(rates) <- sim_conditions
      typeI_step1C <- colSums(rates)
      
      ## Beta 1s
      # Bias
      onestepA_Beta_1i <- rbind(
        PRB_c1 = bias(Beta_1s$class1, parameter = 0.25, type = 'relative', percent = T),
        AB_c1 = bias(Beta_1s$class1, parameter = 0.25, abs = T),
        RSE_c1 = RSE(SE_1s$class1, Beta_1s$class1),
        PRB_c2 = bias(Beta_1s$class2, parameter = -0.25, type = 'relative', percent = T),
        AB_c2 = bias(Beta_1s$class2, parameter = -0.25, abs = T),
        RSE_c2 = RSE(SE_1s$class2, Beta_1s$class2) 
      ) 
      colnames(onestepA_Beta_1i) <- sim_conditions
      step1C_1s <- onestepA_Beta_1i
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1s$class1[i,j] - 1.96*SE_1s$class1[i, j]
          CIs$beta_c1[i, 2] <- Beta_1s$class1[i,j] + 1.96*SE_1s$class1[i, j]
          CIs$beta_c2[i, 1] <- Beta_1s$class2[i,j] - 1.96*SE_1s$class2[i, j]
          CIs$beta_c2[i, 2] <- Beta_1s$class2[i,j] + 1.96*SE_1s$class2[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(0.25))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(-0.25))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      CR <- list("class1" = SEP_c1, "class2" = SEP_c2)
      step1C_1s_cr <- CR
      
      # MSE
      MSE <- rbind(class1 = SimDesign::RMSE(Beta_1s$class1, 0.25, MSE = T),
                   class2 = SimDesign::RMSE(Beta_1s$class2, -0.25, MSE = T))
      colnames(MSE) <- sim_conditions
      step1C_1s_mse <- MSE
    }
    
    if (model== "Step2BC"){

      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 8]
        
      } 
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step2BC_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      step2bc_cvg <-CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      step2BC_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      step2bc_mse <-MSE
    }
    
    if (model== "Step3_A"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 8]
        
      }
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step3A <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      step3a_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      step3A_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      step3A_mse <- MSE
     
    }
    
    if (model== "Step3_B"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 8]
        
      } 
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step3B <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      step3b_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      step3B_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      step3B_mse <- MSE
      
    }
    
    if (model== "Step3_C"){
      ###### Create empty matrix 
      ##### set number of simulation and conditions ####
      nsim <- 100
      conditions <- 8
      
      ####  parameters  ####
      Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      
      Beta_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient of latent intercept
      Beta_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      Beta_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The coefficient for latent slope
      
      theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                       slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      
      Wald_1s <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      
      #### standard error ####
      SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
      SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
      SE_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X2)
      SE_0i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept of latent intercept
      SE_1i <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The slope of latent intercept
      SE_0s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The intercept for latent slope
      SE_1s <- list( class1 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                     class2 = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # The  slope of latent intercept
      SE_theta_3 <- list(); timepoints <- 3
      for (n in 1:timepoints){
        for (indicator in 1:3) {
          SE_theta_3[[paste0("item", n, ".", indicator)]] <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
        } 
      }  # The residual variance for iterms
      
      SE_Sigma_1 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),  # (residual variance for latent intercept)
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))), # (covariance for growth factors)
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))) # (residual variance for latent slope)
      
      SE_Sigma_2 <- list( int = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          int_slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))),
                          slo = matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))))
      for (i in 1:length(results)){
        
        ### Parameters  
        # The intercept of latent class
        Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        Beta_2c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 8]
        
        ### Standard error
        # The intercept of latent class
        SE_0c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 6]
        
        # The slope of latent class (covariates X1)
        SE_1c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 7]
        
        # The slope of latent class (covariates X2)
        SE_2c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 8]
        
      } 
      # Bias
      onestepB_Beta_c <- rbind(
        PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
        AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
        RSE_beta_1c = RSE(SE_1c, Beta_1c),
        PRB_beta_2c = bias(Beta_2c, parameter = 0.75, type = 'relative', percent = T),
        AB_beta_2c = bias(Beta_2c, parameter = 0.75, abs = T),
        RSE_beta_2c = RSE(SE_2c, Beta_2c) 
      ) 
      colnames(onestepB_Beta_c ) <- sim_conditions
      step3C_bias <- onestepB_Beta_c
      
      # CR
      CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)),
                  beta_c2 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
      
      for (j in 1:length(sim_conditions)) {
        CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))),
                    beta_c2 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%"))))
        for(i in 1:nsim){
          CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
          CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
          CIs$beta_c2[i, 1] <- Beta_2c[i,j] - 1.96*SE_2c[i, j]
          CIs$beta_c2[i, 2] <- Beta_2c[i,j] + 1.96*SE_2c[i, j]    
        } # calculate confidence interval for each simulation of conditions
        CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
        CR$beta_c2[j] <- ECR(CIs$beta_c2, parameter = c(0.75))
      } # calculate the coverage rates for each condition
      SEP_c1 <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
      SEP_c2 <- t(BinomCI(CR$beta_c2*100, n=100, conf.level = 0.95, method="wilson"))
      colnames(SEP_c1) <- sim_conditions ; colnames(SEP_c2) <- sim_conditions
      step3c_cvg <- CR
      CR <- list("beta0" = SEP_c1, "beta1" = SEP_c2)
      step3C_cr <- CR
      
      # MSE
      MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T),
                   Beta_2c = SimDesign::RMSE(Beta_2c, 0.75, MSE = T))
      colnames(MSE) <- sim_conditions 
      step3C_mse <- MSE
      
    }
  }

##### compute results #####

  # beta_1c and beta_2c
  #### percentage relative bias ####
table_prb_c1 <- rbind("One-step A" = OnestepA_bias[1, ],
                        "One-step B" = OnestepB_bias[1, ],
                        "One-step C" = OnestepC_bias[1, ],
                        "Two-step A" = step2A1_bias[1, ],
                        "Two-step B" = step2A2B_bias[1, ],
                        "Two-step C" = step2BC_bias[1, ],
                        "Three-step A" = step3A[1, ],
                        "Three-step B" = step3B[1, ],
                        "Three-step C" = step3C_bias[1, ])
  
  table_prb_c2 <- rbind("One-step A" = OnestepA_bias[4, ],
                        "One-step B" = OnestepB_bias[4, ],
                        "One-step C" = OnestepC_bias[4, ],
                        "Two-step A" = step2A1_bias[4, ],
                        "Two-step B" = step2A2B_bias[4, ],
                        "Two-step C" = step2BC_bias[4, ],
                        "Three-step A" = step3A[4, ],
                        "Three-step B" = step3B[4, ],
                        "Three-step C" = step3C_bias[4, ])
  
  
  #### absolute bias ####
  table_ab_c1 <- rbind("One-step A" = OnestepA_bias[2, ],
                       "One-step B" = OnestepB_bias[2, ],
                       "One-step C" = OnestepC_bias[2, ],
                       "Two-step A" = step2A1_bias[2, ],
                       "Two-step B" = step2A2B_bias[2, ],
                       "Two-step C" = step2BC_bias[2, ],
                       "Three-step A" = step3A[2, ],
                       "Three-step B" = step3B[2, ],
                       "Three-step C" = step3C_bias[2, ])
  
  table_ab_c2 <- rbind("One-step A" = OnestepA_bias[5, ],
                       "One-step B" = OnestepB_bias[5, ],
                       "One-step C" = OnestepC_bias[5, ],
                       "Two-step A" = step2A1_bias[5, ],
                       "Two-step B" = step2A2B_bias[5, ],
                       "Two-step C" = step2BC_bias[5, ],
                       "Three-step A" = step3A[5, ],
                       "Three-step B" = step3B[5, ],
                       "Three-step C" = step3C_bias[5, ])
  
  #### standard error ####
  table_se_c1 <- rbind("One-step A" = OnestepA_bias[3, ],
                       "One-step B" = OnestepB_bias[3, ],
                       "One-step C" = OnestepC_bias[3, ],
                       "Two-step A" = step2A1_bias[3, ],
                       "Two-step B" = step2A2B_bias[3, ],
                       "Two-step C" = step2BC_bias[3, ],
                       "Three-step A" = step3A[3, ],
                       "Three-step B" = step3B[3, ],
                       "Three-step C" = step3C_bias[3, ]) %>% round(digits = 2)
  
  
  table_se_c2 <- rbind("One-step A" = OnestepA_bias[6, ],
                       "One-step B" = OnestepB_bias[6, ],
                       "One-step C" = OnestepC_bias[6, ],
                       "Two-step A" = step2A1_bias[6, ],
                       "Two-step B" = step2A2B_bias[6, ],
                       "Two-step C" = step2BC_bias[6, ],
                       "Three-step A" = step3A[6, ],
                       "Three-step B" = step3B[6, ],
                       "Three-step C" = step3C_bias[6, ]) %>% round(digits = 2)
  
  
  
  
  
  #### cr ####
  table_cr_c1lw <- rbind("One-step A" = OnestepA_cr$class1[2,],
                       "One-step B" = OnestepB_cr$beta0[2,],
                       "One-step C" = OnestepC_cr$beta0[2,],
                       "Two-step A" = SteptwoA1_cr$beta0[2,],
                       "Two-step B" = step2A2B_cr$beta0[2,],
                       "Two-step C" = step2BC_cr$beta0[2,],
                       "Three-step A" = step3A_cr$beta0[2,],
                       "Three-step B" = step3B_cr$beta0[2,],
                       "Three-step C" = step3C_cr$beta0[2,])
  table_cr_c1up <- rbind("One-step A" = OnestepA_cr$class1[3,],
                        "One-step B" = OnestepB_cr$beta0[3,],
                        "One-step C" = OnestepC_cr$beta0[3,],
                        "Two-step A" = SteptwoA1_cr$beta0[3,],
                        "Two-step B" = step2A2B_cr$beta0[3,],
                        "Two-step C" = step2BC_cr$beta0[3,],
                        "Three-step A" = step3A_cr$beta0[3,],
                        "Three-step B" = step3B_cr$beta0[3,],
                        "Three-step C" = step3C_cr$beta0[3,])
  
  table_cr_c2lw <- rbind("One-step A" = OnestepA_cr$class2[2,],
                       "One-step B" = OnestepB_cr$beta1[2,],
                       "One-step C" = OnestepC_cr$beta1[2,],
                       "Two-step A" = SteptwoA1_cr$beta1[2,],
                       "Two-step B" = step2A2B_cr$beta1[2,],
                       "Two-step C" = step2BC_cr$beta1[2,],
                       "Three-step A" = step3A_cr$beta1[2,],
                       "Three-step B" = step3B_cr$beta1[2,],
                       "Three-step C" = step3C_cr$beta1[2,])
  
  table_cr_c2up <- rbind("One-step A" = OnestepA_cr$class2[3,],
                         "One-step B" = OnestepB_cr$beta1[3,],
                         "One-step C" = OnestepC_cr$beta1[3,],
                         "Two-step A" = SteptwoA1_cr$beta1[3,],
                         "Two-step B" = step2A2B_cr$beta1[3,],
                         "Two-step C" = step2BC_cr$beta1[3,],
                         "Three-step A" = step3A_cr$beta1[3,],
                         "Three-step B" = step3B_cr$beta1[3,],
                         "Three-step C" = step3C_cr$beta1[3,])
  
  table_cvg_c1 <- rbind("One-step A" = onestepa_cvg$beta_c1,
                        "One-step B" = onestepb_cvg$beta_c1,
                        "One-step C" = onestepc_cvg$beta_c1,
                        "Two-step A" = steptwoa1_cvg$beta_c1,
                        "Two-step B" = step2a2b_cvg$beta_c1,
                        "Two-step C" = step2bc_cvg$beta_c1,
                        "Three-step A" = step3a_cvg$beta_c1,
                        "Three-step B" = step3b_cvg$beta_c1,
                        "Three-step C" = step3c_cvg$beta_c1)
  
  table_cvg_c2 <- rbind("One-step A" = onestepa_cvg$beta_c2,
                        "One-step B" = onestepb_cvg$beta_c2,
                        "One-step C" = onestepc_cvg$beta_c2,
                        "Two-step A" = steptwoa1_cvg$beta_c2,
                        "Two-step B" = step2a2b_cvg$beta_c2,
                        "Two-step C" = step2bc_cvg$beta_c2,
                        "Three-step A" = step3a_cvg$beta_c2,
                        "Three-step B" = step3b_cvg$beta_c2,
                        "Three-step C" = step3c_cvg$beta_c2)
  #### mse ####
  table_MSE_c1 <- rbind("One-step A" = OnestepA_mse[1,],
                        "One-step B" = OnestepB_mse[1,],
                        "One-step C" = OnestepC_mse[1,],
                        "Two-step A" = SteptwoA1_mse[1,],
                        "Two-step B" = step2A2B_mse[1,],
                        "Two-step C" = step2bc_mse[1,],
                        "Three-step A" = step3A_mse[1,],
                        "Three-step B" = step3B_mse[1,],
                        "Three-step C" = step3C_mse[1,])
  
  table_MSE_c2 <- rbind("One-step A" = OnestepA_mse[2,],
                        "One-step B" = OnestepB_mse[2,],
                        "One-step C" = OnestepC_mse[2,],
                        "Two-step A" = SteptwoA1_mse[2,],
                        "Two-step B" = step2A2B_mse[2,],
                        "Two-step C" = step2bc_mse[2,],
                        "Three-step A" = step3A_mse[2,],
                        "Three-step B" = step3B_mse[2,],
                        "Three-step C" = step3C_mse[2,])
  
  
  ##### beta_1i ##### 
  #### bias and SE ####
  prb_1i_c1 <- rbind("One-step B" = OnestepB_1i[1,],
                     "One-step C" = OnestepC_1i[1,],
                     "Step-1 B" = step1B_1i[1,],
                     "Step-1 C" = step1C_1i[1,])
  ab_1i_c1 <- rbind("One-step B" = OnestepB_1i[2,],
                    "One-step C" = OnestepC_1i[2,],
                    "Step-1 B" = step1B_1i[2,],
                    "Step-1 C" = step1C_1i[2,])
  se_1i_c1 <- rbind("One-step B" = OnestepB_1i[3,],
                    "One-step C" = OnestepC_1i[3,],
                    "Step-1 B" = step1B_1i[3,],
                    "Step-1 C" = step1C_1i[3,])
  
  prb_1i_c2 <- rbind("One-step B" = OnestepB_1i[4,],
                     "One-step C" = OnestepC_1i[4,],
                     "Step-1 B" = step1B_1i[4,],
                     "Step-1 C" = step1C_1i[4,])
  ab_1i_c2 <- rbind("One-step B" = OnestepB_1i[5,],
                    "One-step C" = OnestepC_1i[5,],
                    "Step-1 B" = step1B_1i[5,],
                    "Step-1 C" = step1C_1i[5,])
  se_1i_c2 <- rbind("One-step B" = OnestepB_1i[6,],
                    "One-step C" = OnestepC_1i[6,],
                    "Step-1 B" = step1B_1i[6,],
                    "Step-1 C" = step1C_1i[6,])
  
  #### CR ####
  
cr_c1_lw <- rbind("One-step B" = onestepB_1i_cr$class1[2,],
                 "One-step C" = onestepC_1i_cr$class1[2,],
                 "Step-1 B" = steponeB_1i_cr$beta0[2,],
                 "Step-1 C" = steponeC_1i_cr$class1[2,])
  
cr_c1_up <- rbind("One-step B" = onestepB_1i_cr$class1[3,],
                    "One-step C" = onestepC_1i_cr$class1[3,],
                    "Step-1 B" = steponeB_1i_cr$beta0[3,],
                    "Step-1 C" = steponeC_1i_cr$class1[3,])  
  
cr_c2_lw <- rbind("One-step B" = onestepB_1i_cr$class2[2,],
                 "One-step C" = onestepC_1i_cr$class2[2,],
                 "Step-1 B" = steponeB_1i_cr$beta1[2,],
                 "Step-1 C" = steponeC_1i_cr$class2[2,])
cr_c2_up <- rbind("One-step B" = onestepB_1i_cr$class2[3,],
                  "One-step C" = onestepC_1i_cr$class2[3,],
                  "Step-1 B" = steponeB_1i_cr$beta1[3,],
                  "Step-1 C" = steponeC_1i_cr$class2[3,])
  
  #### MSE ####
  mse_c1 <- rbind("One-step B" = onestepB_1i_mse[1,],
                  "One-step C" = onestepC_1i_mse[1,],
                  "Step-1 B" = steponeB_1i_mse[1,],
                  "Step-1 C" = steponeC_1i_mse[1,])
  
  mse_c2 <- rbind("One-step B" = onestepB_1i_mse[2,],
                  "One-step C" = onestepC_1i_mse[2,],
                  "Step-1 B" = steponeB_1i_mse[2,],
                  "Step-1 C" = steponeC_1i_mse[2,])
  

  
  
  
  ##### beta_1s #####
  ### bias and SE ####
  prb_1s_c1 <- rbind("One-step C" = OnestepC_1s[1,],
                     "Step-1 C" = step1C_1s[1,])
  ab_1s_c1 <- rbind("One-step C" = OnestepC_1s[2,],
                    "Step-1 C" = step1C_1s[2,])
  se_1s_c1 <- rbind("One-step C" = OnestepC_1s[3,],
                    "Step-1 C" = step1C_1s[3,])
  
  prb_1s_c2 <- rbind("One-step C" = OnestepC_1s[4,],
                     "Step-1 C" = step1C_1s[4,])
  ab_1s_c2 <- rbind("One-step C" = OnestepC_1s[5,],
                    "Step-1 C" = step1C_1s[5,])
  se_1s_c2 <- rbind("One-step C" = OnestepC_1s[6,],
                    "Step-1 C" = step1C_1s[6,])
  ### CR ####
  cr_1s_c1lw <- rbind("One-step C" = onestepC_1s_cr$class1[2,],
                 "Step-1 C" = step1C_1s_cr$class1[2,])
  cr_1s_c1up <- rbind("One-step C" = onestepC_1s_cr$class1[3,],
                      "Step-1 C" = step1C_1s_cr$class1[3,])
  cr_1s_c2lw <- rbind("One-step C" = onestepC_1s_cr$class2[2,],
                 "Step-1 C" = step1C_1s_cr$class2[2,])
  cr_1s_c2up <- rbind("One-step C" = onestepC_1s_cr$class2[3,],
                      "Step-1 C" = step1C_1s_cr$class2[3,])
  
  #### MSE ####
  mse_1s_c1 <- rbind("One-step C" = onestepC_1s_mse[1,],
                  "Step-1 C" = step1C_1s_mse[1,])
  
  mse_1s_c2 <- rbind("One-step C" = onestepC_1s_mse[2,],
                  "Step-1 C" = step1C_1s_mse[2,])
  
  
  
  
  
  
  
  
  
  
  
#### storage the result ####
results_3t <- list(table_ab_c1=table_ab_c1, table_ab_c2=table_ab_c2,
                   table_cr_c2up=table_cr_c2up, table_cr_c1up=table_cr_c1up,
                   table_cr_c2lw=table_cr_c2lw, table_cr_c1lw=table_cr_c1lw,
                   table_cvg_c1=table_cvg_c1, table_cvg_c2=table_cvg_c2,
                   table_MSE_c1=table_MSE_c1, table_MSE_c2=table_MSE_c2, 
                   table_prb_c1=table_prb_c1, table_prb_c2=table_prb_c2,
                   table_se_c1=table_se_c1, table_se_c2=table_se_c2,
                   
                   prb_1i_c1=prb_1i_c1, prb_1i_c2=prb_1i_c2,
                   se_1i_c1=se_1i_c1, se_1i_c2=se_1i_c2,
                   cr_c1_lw=cr_c1_lw, cr_c2_up=cr_c2_up,
                   cr_c1_up=cr_c1_up, cr_c2_lw=cr_c2_lw,
                   ab_1i_c1=ab_1i_c1, ab_1i_c2=ab_1i_c2,
                   mse_c1=mse_c1, mse_c2=mse_c2,
                   typeI = rbind(typeI_onestepC,typeI_step1C),
                   
                   prb_1s_c1 = prb_1s_c1, prb_1s_c2 = prb_1s_c2,
                   ab_1s_c1 = ab_1s_c1, ab_1s_c2 = ab_1s_c2,
                   se_1s_c1 = se_1s_c1, se_1s_c2 = se_1s_c2,
                   cr_1s_c1lw = cr_1s_c1lw, cr_1s_c2up = cr_1s_c2up,
                   cr_1s_c1up = cr_1s_c1up, cr_1s_c2lw = cr_1s_c2lw,
                   mse_1s_c1 = mse_1s_c1, mse_1s_c2 = mse_1s_c2
                   )

results_6t <- list(table_ab_c1=table_ab_c1, table_ab_c2=table_ab_c2,
                   table_cr_c2up=table_cr_c2up, table_cr_c1up=table_cr_c1up,
                   table_cr_c2lw=table_cr_c2lw, table_cr_c1lw=table_cr_c1lw,
                   table_cvg_c1=table_cvg_c1, table_cvg_c2=table_cvg_c2,
                   table_MSE_c1=table_MSE_c1, table_MSE_c2=table_MSE_c2, 
                   table_prb_c1=table_prb_c1, table_prb_c2=table_prb_c2,
                   
                   table_se_c1=table_se_c1, table_se_c2=table_se_c2,
                   prb_1i_c1=prb_1i_c1, prb_1i_c2=prb_1i_c2,
                   se_1i_c1=se_1i_c1, se_1i_c2=se_1i_c2,
                   cr_c1_lw=cr_c1_lw, cr_c2_up=cr_c2_up,
                   cr_c1_up=cr_c1_up, cr_c2_lw=cr_c2_lw,
                   ab_1i_c1=ab_1i_c1, ab_1i_c2=ab_1i_c2,
                   mse_c1=mse_c1, mse_c2=mse_c2,
                   typeI = rbind(typeI_onestepC,typeI_step1C),
                   
                   prb_1s_c1 = prb_1s_c1, prb_1s_c2 = prb_1s_c2,
                   ab_1s_c1 = ab_1s_c1, ab_1s_c2 = ab_1s_c2,
                   se_1s_c1 = se_1s_c1, se_1s_c2 = se_1s_c2,
                   cr_1s_c1lw = cr_1s_c1lw, cr_1s_c2up = cr_1s_c2up,
                   cr_1s_c1up = cr_1s_c1up, cr_1s_c2lw = cr_1s_c2lw,
                   mse_1s_c1 = mse_1s_c1, mse_1s_c2 = mse_1s_c2
                  )
study2_3t <-results_3t
study2_6t <-results_6t

study3_3t <-results_3t
study3_6t <-results_6t

#### average time points ####
mean_results <- results_6t
for (i in 1:length(results_3t)) {
  mean_results[[i]] <- (results_3t[[i]]+results_6t[[i]])/2
}

####################################### beta_xi and beta_x2
#### table for Bias, CVG and SE/SD #####
merge_table <- function(Bias, Coverage_rate_lw, Coverage_rate_up,SE) {
  
  format_cell <- function(value1, value2, value3) {
    ifelse( value3 >= 0.95,
           value2 <- sprintf("%.2f (%.2f-%.2f)", value1, value2, value3),
           value2 <- sprintf("%.2f (\\textbf{%.2f-%.2f})", value1, value2, value3))
  }
  
  
  formatted_data <- Bias
  for(i in 1:ncol(Bias)){
    for(j in 1:nrow(Bias)){
      formatted_data[j,i] <- format_cell(Bias[j,i], Coverage_rate_lw[j,i], Coverage_rate_up[j,i])
    }
  }
 
  # 将两列表格交叉合并
  merged_data <- matrix(NA,nrow = nrow(formatted_data), ncol = ncol(Bias) + ncol(SE),
                        dimnames=list(c(rownames(Bias)),c(colnames(SE),colnames(SE))))
  for (i in 1:(ncol(formatted_data) + ncol(SE))) {
    if (i %% 2 != 0) {
      merged_data[, i] <- formatted_data[, (i + 1) %/% 2]
    } else {
      merged_data[, i] <- SE[, i %/% 2]
    }
  }
  
  # Create a table using kable 
  Table <-  kable(merged_data, "latex", booktabs = T, escape = F,
                  caption = "The Percentage Relative Bias (Coverage) Values Over all Replications for the Latent Class Parameters $\\beta_{1x_1}$ and $\\beta_{1x_2}$, in Eight Simulation Conditions \\label{tab:beta_c bias}", 
                  col.names = c(rep(c("$Bias$($CVG$)", "$SE/SD$"), 8)), align = "c") %>% 
    add_header_above(c(" "=1, "$N=1000$" = 2, "$N=500$" = 2,"$N=1000$" = 2, "$N=500$" = 2, "$N=1000$" = 2, "$N=500$" = 2, "$N=1000$" = 2, "$N=500$" = 2),escape = F) %>%
    add_header_above(c(" "=1, "High entropy" = 4, "Moderate entropy" = 4, "High entropy" = 4, "Moderate entropy" = 4 )) %>%
    add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 8, "Mixing ratio = 0.30/0.70" = 8)) %>%
    pack_rows(index=c("$\\beta_{1x_1}$ = -0.5" = 9, "$\\beta_{1x_2}$ = 0.75" = 9),escape = F) %>% 
    footnote(general = "Two step (Three step) A = Two-step (Three-step) estimator ignore the DE at step 1 model. Two step (Three step) B = Two-step (Three-step) estimator only model the DE on the latent intercept at step 1 model. Two step (Three step) C = Two-step (Three-step) estimator model the DE on both growth factors at step 1 model. " , threeparttable = T) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 7) %>%
    landscape()
  
  print(Table)
}

ab <- rbind(mean_results$table_ab_c1, mean_results$table_ab_c2)
cr_lw <- rbind(mean_results$table_cr_c1lw, mean_results$table_cr_c2lw)
cr_up <- rbind(mean_results$table_cr_c1up, mean_results$table_cr_c2up)
se <- rbind(mean_results$table_se_c1, mean_results$table_se_c2)
merge_table(ab, cr_lw, cr_up, se)
####



#### table for Bias and CVG ####
format_table <- function(Bias, Coverage_rate_lw, Coverage_rate_up) {
  
  
  format_cell <- function(value1, value2, value3) {
    ifelse( value3 >= 0.95 & value2 <= 0.95,
           value2 <- sprintf("%.2f (%.2f-%.2f)", value1, value2, value3),
           value2 <- sprintf("%.2f (\\textbf{%.2f-%.2f})", value1, value2, value3))
  }
  
  formatted_data <- Bias  
  for(i in 1:ncol(Bias)){
    for(j in 1:nrow(Bias)){
      formatted_data[j,i] <- format_cell(Bias[j,i], Coverage_rate_lw[j,i], Coverage_rate_up[j,i])
    }
  }
  
  
  # Create a table using kable 
  Table <-  kable(formatted_data, "latex", booktabs = T, escape = F,
                  caption = "The Absolute Bias ($95\\%$ Confidence Interval of Coverage Rates) Values Over all Replications 
                  for the $\\beta_{x_1}$ and $\\beta_{x_2}$, in 6 Time Points 
                  Conditions for Study 2 \\label{tab:beta_x1x2 study2}", 
                  col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
    add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
    add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
    pack_rows(index=c("$\\beta_{x_1}$ = -0.5" = 9, "$\\beta_{x_2}$ = 0.75" = 9),escape = F) %>% 
    footnote(general = "Two-step (Three-step) A = Two-step (Three-step) estimator ignores the DE at step 1 model. 
             Two-step (Three-step) B = Two-step (Three-step) estimator only models the DE on the latent intercept 
             at step 1 model. Two-step (Three-step) C = Two-step (Three-step) estimator model the DE on both growth 
             factors at step 1 model. $N$ = total sample size. The bold numbers reflect the coverage rates are lower 
             than $95\\%$ in these conditions. " , threeparttable = T,escape = F) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 7)# %>%
    #cell_spec(formatted_data[2,1],bold = T) 
    #%>% landscape()
   
  print(Table)
}
ab <- rbind(mean_results$table_ab_c1, mean_results$table_ab_c2)
cr_lw <- rbind(mean_results$table_cr_c1lw, mean_results$table_cr_c2lw)
cr_up <- rbind(mean_results$table_cr_c1up, mean_results$table_cr_c2up)
se <- rbind(mean_results$table_se_c1, mean_results$table_se_c2)
formatted_data <- format_table(ab, cr_lw,cr_up)



#### table for MSE ####
MSE_1c <- rbind("beta c1" = mean_results$table_MSE_c1,
                "beta c2" = mean_results$table_MSE_c2)
kable(MSE_1c, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Mean Square Error Values Over all Replications for the $\\beta_{x_1}$ and $\\beta_{x_2}$ in 3 Time Points Conditions \\label{tab:beta_x1x2 mse 3T} ", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\beta_{x_1} = -0.5$" = 9, "$\\beta_{x_2} = 0.75$" = 9), escape =F) %>% 
  footnote(general = "One-step A, One-step B, and One-step C are the one-step estimators ignoring the direct effects (specifications A), specifying the direct effect on the latent intercept (specification B), and on both growth factors (specification C), respectively. Two-step A, Two-step B, and Two-step C are the proposed two-step estimators with specifications A, B, and C, respectively. Three-step A, Three-step B, and Three-step C are the three-step estimators with specifications A, B, and C, respectively. $N$ is the total sample size. The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval. " , threeparttable = T) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)



#### table for SE/SD ####
SE_1c <- rbind("beta c1" = mean_results$table_se_c1,
              "beta c2" = mean_results$table_se_c2)
kable(SE_1c, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Standard Error Ratio Values Over all Replications for the $\\beta_{x_1}$ and $\\beta_{x_2}$ in 3 Time Points Conditions \\label{tab:beta_x1x2 se 3T}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\beta_{x_1} = -0.5$" = 9, "$\\beta_{x_2} = 0.75$" = 9), escape =F) %>% 
  footnote(general = "One-step A, One-step B, and One-step C are the one-step estimators ignoring the direct effects (specifications A), specifying the direct effect on the latent intercept (specification B), and on both growth factors (specification C), respectively. Two-step A, Two-step B, and Two-step C are the proposed two-step estimators with specifications A, B, and C, respectively. Three-step A, Three-step B, and Three-step C are the three-step estimators with specifications A, B, and C, respectively. $N$ is the total sample size. The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval. " , threeparttable = T) %>%
  kable_styling(latex_options = c("hold_position"), font_size =7)


####################################### gamma_1 ######
mean_results <- results_3t

#### table for Bias, CVG and SE/SD #### 
Bias_1i <- rbind(mean_results$ab_1i_c1, mean_results$ab_1i_c2)
cr_lw <- rbind(mean_results$cr_c1_lw, mean_results$cr_c2_lw)
cr_up <- rbind(mean_results$cr_c1_up, mean_results$cr_c2_up)
format_table_1i <- function(Bias, Coverage_rate_lw, Coverage_rate_up) {
  
  format_cell <- function(value1, value2, value3) {
    ifelse( value3 >= 0.95 & value2 <= 0.95,
            value2 <- sprintf("%.2f (%.2f-%.2f)", value1, value2, value3),
            value2 <- sprintf("%.2f (\\textbf{%.2f-%.2f})", value1, value2, value3))
  }
  
  formatted_data <- Bias  
  for(i in 1:ncol(Bias)){
    for(j in 1:nrow(Bias)){
      formatted_data[j,i] <- format_cell(Bias[j,i], Coverage_rate_lw[j,i], Coverage_rate_up[j,i])
    }
  }
  
  # Create a table using kable 
  Table <-  kable(formatted_data, "latex", booktabs = T, escape = F, 
                  caption = "The Absolute Bias (95\\% Confidence Interval of Coverage Rates) 
                  Values Over all Replications for the Latent intercept Parameter $\\gamma_{0}$ 
                  in 6 Time Points Conditions \\label{tab:gamma_0 study2}", 
                  col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
    add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
    add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
    pack_rows(index=c("$\\gamma_0^{(1)}=0.50$" = 4, "$\\gamma_0^{(2)}=-0.50$" = 4),escape = F) %>% 
    footnote(general = "One-step B is the one-step estimator that only models the direct effect (DE) 
             on the latent intercept (specification B). One-step C is the one-step estimator that models 
             the DE on both growth factors (specification C). Step-1 B is the step-1 model of the proposed 
             two-step estimator and the three-step estimator with specification B.  Step-1 C is the step-1 
             model of the proposed two-step estimator and the three-step estimator with specification C. 
             $\\gamma_0^{(1)}$ is the regression coefficients of latent intercept in Class 1. $\\gamma_0^{(2)}$ 
             is the regression coefficients of latent intercept in Class 2. $N$ is the total sample size. 
             The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval.", 
             threeparttable = T,escape = F ) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 8)
  
  print(Table)
}
format_table_1i(Bias_1i, cr_lw, cr_up)

#### table for MSE ####

MSE_1c3 <- rbind(study2_3t$table_MSE_c1, study2_3t$table_MSE_c2)
MSE_1c6 <- rbind(study2_6t$table_MSE_c1, study2_6t$table_MSE_c2)

MSE_1c3 <- rbind(study3_3t$table_MSE_c1, study3_3t$table_MSE_c2)
MSE_1c6 <- rbind(study3_6t$table_MSE_c1, study3_6t$table_MSE_c2)

MSE_1c <- cbind(MSE_1c3,MSE_1c6)
MSE_long <- rbind(MSE_1c[, 1:8], MSE_1c[, 9:16])

# wide table (data: MSE_1c)
kable(MSE_1c, "latex", digits = 3, booktabs = T, escape = F,
      caption = "The Mean Square Error Values Over all Replications for the $\\beta_{x_1}$, and $\\beta_{x_2}$ in Eight Simulation Conditions in study 3 \\label{tab:beta_i mse}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,  
                     "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4, 
                     "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
  pack_rows(index=c("$\\beta_{x_1}$" = 9, "$\\beta_{x_2}$" = 9),escape =F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)

# long table (data: MSE_long)
kable(MSE_long, "latex", digits = 3, booktabs = T, escape = F,
      caption = "The Mean Square Error Values Over all Replications for the $\\beta_{x_1}$, and $\\beta_{x_2}$ in Eight Simulation Conditions in study 3 \\label{tab:beta_i mse}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  #add_header_above(c(" "=1, "Three time points" = 4, "Six time points" = 4)) %>%
  pack_rows(index = c("Three time points" = 18, "Six time points" = 18)) %>%
  pack_rows(index = c("$\\beta_{x_1}$" = 9, "$\\beta_{x_2}$" = 9,"$\\beta_{x_1}$" = 9, "$\\beta_{x_2}$" = 9),escape =F) %>% 
  footnote(general = "One-step B is the one-step estimator that only models the direct effect (DE) 
             on the latent intercept (specification B). One-step C is the one-step estimator that models 
             the DE on both growth factors (specification C). Step-1 B is the step-1 model of the proposed 
             two-step estimator and the three-step estimator with specification B.  Step-1 C is the step-1 
             model of the proposed two-step estimator and the three-step estimator with specification C. 
             $\\gamma_0^{(1)}$ is the regression coefficients of latent intercept in Class 1. $\\gamma_0^{(2)}$ 
             is the regression coefficients of latent intercept in Class 2. $N$ is the total sample size. 
             The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval.", 
           escape = F ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)



#### single factor tables ####

se_13 <- rbind(study2_3t$table_se_c1,study2_3t$table_se_c2)
se_16 <- rbind(study2_6t$table_se_c1,study2_6t$table_se_c2)

se_13 <- rbind(study3_3t$table_se_c1,study3_3t$table_se_c2)
se_16 <- rbind(study3_6t$table_se_c1,study3_6t$table_se_c2)
se_1c <- cbind(se_13, se_16)
se_1c_long <- rbind(se_13, se_16)

kable(se_1c, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Standard Error Ratio Values Over all Replications for the Latent Class Parameters $\\gamma_0$ in Study 3 \\label{tab:beta_i se_sig_stu3}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,  "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
  pack_rows(index=c("$\\beta_{x_1}}$" = 9, "$\\beta_{x_2}}$" = 9),escape =F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)

ab_1s <- rbind("beta c1" = mean_results$ab_1i_c1,
               "beta c2" = mean_results$ab_1i_c2)
kable(ab_1s, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Absolute Bias Values Over all Replications for the Latent Class Parameters $\\beta_{X1i}$ in Eight Simulation Conditions \\label{tab:beta_i se}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("Class 1" = 4, "Class 2" = 4)) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)


#### table for type I rates #### 
rownames(mean_results$typeI) <- c("One step C", "Step one C")
kable(mean_results$typeI*0.01, "latex",digits = 2, booktabs = T, escape = F, 
      caption = "The Type I Error Rate Values Over all Replications for the Latent Slope Regression Coefficients $\\bm{\\gamma}_1$ in 6 Time Points Conditions for Study 2  \\label{tab:beta_s type I}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  footnote(general = "One-step C is the one-step estimator model of the direct effects (DEs) on both growth factors (specification C). Step-1 C is the step-1 model of the two-step and three-step estimators with specification C. $N$ is the total sample size.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 8)




#### beta 1s ####
Bias_1s <- rbind(mean_results$ab_1s_c1, mean_results$ab_1s_c2)
cr_1s_lw <- rbind(mean_results$cr_1s_c1lw, mean_results$cr_1s_c2lw)
cr_1s_up <- rbind(mean_results$cr_1s_c1up, mean_results$cr_1s_c2up)
format_table_1i <- function(Bias, Coverage_rate_lw, Coverage_rate_up) {
  
  format_cell <- function(value1, value2, value3) {
    ifelse( value3 >= 0.95 & value2 <= 0.95,
            value2 <- sprintf("%.2f (%.2f-%.2f)", value1, value2, value3),
            value2 <- sprintf("%.2f (\\textbf{%.2f-%.2f})", value1, value2, value3))
  }
  
  formatted_data <- Bias  
  for(i in 1:ncol(Bias)){
    for(j in 1:nrow(Bias)){
      formatted_data[j,i] <- format_cell(Bias[j,i], Coverage_rate_lw[j,i], Coverage_rate_up[j,i])
    }
  }
  
  # Create a table using kable 
  Table <-  kable(formatted_data, "latex", booktabs = T, escape = F, 
                  caption = "The Absolute Bias (95\\% Confidence Interval of Coverage Rates) Values Over all Replications
                  for the Latent intercept Parameter $\\gamma_1$ in 6 Time Points Conditions for Study 3 \\label{tab:gamma_1 bias}", 
                  col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
    add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
    add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
    pack_rows(index=c("Class 1" = 2, "Class 2" = 2)) %>% 
    footnote(general = "One-step C is the one-step estimator that specifies direct effects (DEs) on both growth factors (specification C). 
             Step-1 C is the step-1 model of the two-step and the three-step estimators that specifies the direct effects 
             on growth factors (specification C).", threeparttable = T ) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 8)
  
  print(Table)
}
format_table_1i(Bias_1s,cr_1s_lw, cr_1s_up)

# MSE for gamma0

MSE_1i3 <- rbind(study2_3t$mse_c1, study2_3t$mse_c2)
MSE_1i6 <- rbind(study2_6t$mse_c1, study2_6t$mse_c2)

MSE_1i3 <- rbind(study3_3t$mse_c1, study3_3t$mse_c2)
MSE_1i6 <- rbind(study3_6t$mse_c1, study3_6t$mse_c2)
MSE_1i <- rbind(MSE_1i3, MSE_1i6) # gamma_0
rownames(SE_1i) <- c(rep(c("One-step B", "One-step C", "step-one B", "step-one c"), 4))

MSE_1s3 <- rbind(study2_3t$mse_1s_c1, study2_3t$mse_1s_c2)
MSE_1s6 <- rbind(study2_6t$mse_1s_c1, study2_6t$mse_1s_c2)

MSE_1s3 <- rbind(study3_3t$mse_1s_c1, study3_3t$mse_1s_c2)
MSE_1s6 <- rbind(study3_6t$mse_1s_c1, study3_6t$mse_1s_c2)
MSE_1s <- rbind(MSE_1s3, MSE_1s6) # gamma_1
rownames(SE_1s) <- c(rep(c( "One-step C", "step-one c"), 4))

kable(MSE_1s, "latex", digits = 3, booktabs = T, escape = F,
      caption = "The Mean Square Error Values Over all Replications for the $\\gamma_0$ in Eight Simulation Conditions \\label{tab:beta_i mse}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,  "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
  pack_rows(index=c("$\\gamma_1^{1}=0.25$" = 2, "$\\gamma_1^{2}=-0.25$" = 2),escape =F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) 
           only model the DE on the latent intercept. One step (step one) C = one-step estimator 
           (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)


# SE
SE_1i3 <- rbind(study2_3t$se_1i_c1, study2_3t$se_1i_c2)
SE_1i6 <- rbind(study2_6t$se_1i_c1, study2_6t$se_1i_c2)

SE_1i3 <- rbind(study3_3t$se_1i_c1, study3_3t$se_1i_c2)
SE_1i6 <- rbind(study3_6t$se_1i_c1, study3_6t$se_1i_c2)
SE_1i <- rbind(SE_1i3, SE_1i6) # gamma_0
rownames(SE_1i) <- c(rep(c("One-step B", "One-step C", "step-one B", "step-one c"), 4))

SE_1s3 <- rbind(study2_3t$se_1s_c1, study2_3t$se_1s_c2)
SE_1s6 <- rbind(study2_6t$se_1s_c1, study2_6t$se_1s_c2)

SE_1s3 <- rbind(study3_3t$se_1s_c1, study3_3t$se_1s_c2)
SE_1s6 <- rbind(study3_6t$se_1s_c1, study3_6t$se_1s_c2)
SE_1s <- rbind(SE_1s3, SE_1s6) # gamma_1
rownames(SE_1s) <- c(rep(c( "One-step C", "step-one c"), 4))



se_13 <- rbind(results_3t$se_1i_c1,results_3t$se_1i_c2)
se_16 <- rbind(results_6t$se_1i_c1,results_6t$se_1i_c2)
se_1i <- cbind(se_13, se_16)

kable(se_1i, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Standard Error Ratio Values Over all Replications for the $\\gamma_0$ in Study 3 \\label{tab:beta_i se}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,"High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4,"Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\gamma_0^{(1)}=0.25$" = 4, "$\\gamma_0^{(2)}=-0.25$" = 4),escape = F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)

se_13 <- rbind(results_3t$se_1s_c1,results_3t$se_1s_c2)
se_16 <- rbind(results_6t$se_1s_c1,results_6t$se_1s_c2)
se_1s <- cbind(se_13, se_16)

kable(se_1s, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Standard Error Ratio Values Over all Replications for the $\\gamma_1$ in Study 3 \\label{tab:beta_i se}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,"High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4,"Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\gamma_1^{(1)}=0.25$" = 2, "$\\gamma_1^{(2)}=-0.25$" = 2),escape = F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)



# AB

ab_1s <- rbind("beta c1" = mean_results$ab_1i_c1,
               "beta c2" = mean_results$ab_1i_c2)
kable(ab_1s, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Absolute Bias Values Over all Replications for the Latent Class Parameters $\\beta_{X1i}$ in Eight Simulation Conditions \\label{tab:beta_i se}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\gamma_0^{(1)}=0.25$" = 4, "$\\gamma_0^{(2)}=-0.25$" = 4),escape = F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)


# type I 
rownames(mean_results$typeI) <- c("One step C", "Step one C")
kable(mean_results$typeI*0.01, "latex",digits = 2, booktabs = T, escape = F, 
      caption = "The Type I Error Rate Values Over all Replications for the Latent Slope Parameter $\\beta_{X2s}$ in Eight Simulation Conditions \\label{tab:beta_s type I}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  footnote(general = "One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 8)

#------------------------------------------------------------------------------
###### Graph Part #########
library(ggplot2)
library(tidyverse)
library(tikzDevice)
library(openxlsx)
library(tidyverse)
library(ggsci)
library(latex2exp)

# Turning the rownames into the first column and tranfer to a long version of data  
Long_beta_c_3 <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_3t$table_se_c1), var="Models"),
                                         rownames_to_column(as.data.frame(results_3t$table_se_c2), var="Models")), 
                                    cols = EH1000:UL500,
                                    names_to = "Conditions",
                                    values_to = "SE_ratio"),
                      Par=c(rep(c("beta_1x1"),144/2), rep(c("beta_1x2"),144/2)),
                      Time_points = "Three time points")
Long_beta_c_6 <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_6t$table_se_c1), var="Models"),
                                         rownames_to_column(as.data.frame(results_6t$table_se_c2), var="Models")), 
                                   cols = EH1000:UL500,
                                   names_to = "Conditions",
                                   values_to = "SE_ratio"),
                      Par=c(rep(c("beta_1x1"),144/2), rep(c("beta_1x2"),144/2)),
                      Time_points = "Six time points")
Long_tab <- rbind(Long_beta_c_3,Long_beta_c_6)


# Define the latex version of factor names
f_names <- list(
  "beta_1x1" = TeX(c("$\\beta_{x_1}$")),
  "beta_1x2" = TeX(c("$\\beta_{x_2}$")),
  "Three time points" = "Three time points" ,
  "Six time points" = "Six time points" 
)
f_labeller <- function (variable, value){
  return(f_names[value])
}

# Generate plots
#gapminder %>% filter(year == 2002) %>%
ggplot(Long_tab, aes(y=SE_ratio,x=Models)) +
  geom_boxplot() +
  #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
  #  geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
  facet_grid(vars(Time_points),vars(Par), labeller = f_labeller, scales = "free_x")+
  geom_hline(yintercept = 1,color = "blue")+
  theme_bw() +
  coord_flip()+
  labs(y =TeX(r'($SE/SD)'))+
  theme(axis.text.y = element_text(face=rep(c("plain","bold","plain"),3),size = 10))+ # Study 2 bold specification B
  #theme(axis.text.y = element_text(face=rep(c("bold","plain","plain"),3), size = 12))+ # Study 3 bold specification C
  theme(axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 12))+
  scale_x_discrete(limits=c( "Three-step C", "Three-step B", "Three-step A", 
                             "Two-step C", "Two-step B", "Two-step A", 
                             "One-step C", "One-step B", "One-step A"))



Long_gamma1_3  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_3t$se_1i_c1), var="Models"),
                                           rownames_to_column(as.data.frame(results_3t$se_1i_c2), var="Models")), 
                                    cols = EH1000:UL500,
                                    names_to = "Conditions",
                                    values_to = "SE_ratio" ),
                       Par = c("gamma_1" ),
                       Class = c(rep(c("Class 1"), 32),rep("Class 2", 32)),
                       Time_points = "Three time points")
Long_gamma1_6  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_6t$se_1i_c1) ,var="Models"),
                                           rownames_to_column(as.data.frame(results_6t$se_1i_c2) ,var="Models")), 
                                     cols = EH1000:UL500,
                                     names_to = "Conditions",
                                     values_to = "SE_ratio"),
                        Par=c("gamma_1"),
                        Class = c(rep(c("Class 1"), 32),rep("Class 2", 32)),
                        Time_points = "Six time points")
Long_tab <- rbind(Long_gamma1_3,Long_gamma1_6)

# recode model names
Long_tab <- Long_tab %>%
  mutate(Models = recode(Models,
                         "Step-1 B" = "Step-one B",
                         "Step-1 C" = "Step-one C"))

# Generate plots
#gapminder %>% filter(year == 2002) %>%
ggplot(Long_tab, aes(y=SE_ratio,x=Models)) +
  geom_boxplot() +
  #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
  #  geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
  geom_hline(yintercept = 1,color = "blue")+
  facet_grid(.~Time_points)+
  theme_bw() + coord_flip()+
  labs(y =TeX(r'($SE/SD)'), x = "Models")+
  ylim(0.80,1.30)+
  theme(axis.text.y = element_text(face=rep(c("plain", "bold"),2), size = 16))+ # study 2 bold specification B
#  theme(axis.text.y = element_text(face=rep(c("bold", "plain"),2), size = 16))+ # study 3 bold specification C
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 14))+
  scale_x_discrete(limits=c( "Step-one C", "Step-one B",
                           "One-step C", "One-step B"))
  

# gamma_1
Long_gamma1_3  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_3t$se_1s_c1), var="Models"),
                                           rownames_to_column(as.data.frame(results_3t$se_1s_c2), var="Models")), 
                                     cols = EH1000:UL500,
                                     names_to = "Conditions",
                                     values_to = "SE_ratio" ),
                        Par = c("gamma_1" ),
                        Class = c(rep(c("Class 1"), 32),rep("Class 2", 32)),
                        Time_points = "Three time points")
Long_gamma1_6  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_6t$se_1s_c1) ,var="Models"),
                                           rownames_to_column(as.data.frame(results_6t$se_1s_c2) ,var="Models")), 
                                     cols = EH1000:UL500,
                                     names_to = "Conditions",
                                     values_to = "SE_ratio"),
                        Par=c("gamma_1"),
                        Class = c(rep(c("Class 1"), 32),rep("Class 2", 32)),
                        Time_points = "Six time points")
Long_tab <- rbind(Long_gamma1_3,Long_gamma1_6)
# recode model names
Long_tab <- Long_tab %>%
  mutate(Models = recode(Models,
                         "One step C" = "One-step C",
                         "Step-1 C" = "Step-one C"))

# Generate plots
#gapminder %>% filter(year == 2002) %>%
gamma1 <- ggplot(Long_tab, aes(y=SE_ratio,x=Models)) +
  geom_boxplot() +
  #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
  #  geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
  geom_hline(yintercept = 1,color = "blue")+
  facet_grid(.~Time_points)+
  theme_bw() + coord_flip()+
  labs(y =TeX(r'($SE/SD)'), x = "Models")+
  ylim(0.80,1.20)+
  theme(axis.text.y = element_text(face=rep(c("plain","plain"),3), size = 14))+
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 14))+
  scale_x_discrete(limits=c( "Step-one C", "One-step C"))

library(ggplot2)
library(scales)
require(grid)
library(tikzDevice)
tinytex::tlmgr_install(pkgs = "tikz")
tinytex::install_tinytex()
tinytex::tlmgr_install(pkgs = "pdftex-def")


library(tinytex)
library(filehash)
tinytex::tlmgr_install(pkgs = "tikz")
options(tikzDefaultEngine = 'pdftex')
options(tikzDocumentDeclaration = "\\documentclass[11pt]{report}")

library(tinytex)
tlmgr_search('/tikz.sty')    # search for tikz.sty
tlmgr_install('pdftex-def')         # install the psnfss package
tlmgr_update()               # update everything



tikz(file = "my_output_file.tex", standAlone=F,width = 6, height = 6)

plot <-  ggplot(Long_tab, aes(y=SE_ratio,x=Models)) +
  geom_boxplot() +
  #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
  #  geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
  geom_hline(yintercept = 1,color = "blue")+
  facet_grid(.~Time_points)+
  theme_bw() + coord_flip()+
  labs(y =TeX(r'($SE/SD)'), x = "Models")+
  ylim(0.80,1.20)+
  theme(axis.text.y = element_text(face=rep(c("plain","plain"),3), size = 14))+
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 14))+
  scale_x_discrete(limits=c( "Step-one C", "One-step C"))

        print(plot)
        
        dev.off()

        td <- "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/Figures"
        tf <- file.path(td,'test1.tex')
        oldwd <- getwd()
        setwd(td)
        
        # Minimal plot
        tikz(file = 'SE_study1.tex')
        p <- ggplot(Long_tab, aes(x=SE_ratio,y=Models)) +
          geom_boxplot() +
          #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
          #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
          facet_wrap(vars(Time_points)
                     #, labeller = f_labeller, scales = "free_x"
          )+
          geom_vline(xintercept = 1,color = "blue")+
          theme_bw() +
          xlim(c(0.9,1.3))+
          #coord_flip()+
          labs(x =TeX(r'($SE/SD)'), y= "Models")+
          theme(axis.text.y = element_text(face=rep(c("plain","plain","plain"),3), size = 14))+
          theme(axis.text.x = element_text(size = 12),
                axis.title.y = element_text(size = 16),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 12))+
          scale_y_discrete(limits=c( "Three-step", "Two-step", "One-step"))
        plot(p)
        dev.off()
        
        # View the output
        tools::texi2dvi(plot,pdf=T)
        system(paste(getOption('pdfviewer'),file.path(td,'example1.pdf')))
        setwd(oldwd)
        7
        
        packs <- c("`xcolor.sty'") 
        
        for (i in 1:length(packs)) {     
          tinytex::parse_install(text = paste('! LaTeX Error: File', packs[i],'not found.'))
        }

# Entropy

Entropy_s3_3 <- Entropy_s3
Entropy_s3_6 <- Entropy_s3
Entropy_s2_full <- rbind(Entropy_s2,Entropy_s2_6)
Entropy_s3_full <- rbind(Entropy_s3_3,Entropy_s3_6)
Entropy_low <- rbind(Entropy_s2_full[,c(3,4,7,8)],Entropy_s3_full[,c(3,4,7,8)])
Entropy_high <- rbind(Entropy_s2_full[,c(1,2,5,6)],Entropy_s3_full[,c(1,2,5,6)])
summary(Entropy_low)
summary(Entropy_high)