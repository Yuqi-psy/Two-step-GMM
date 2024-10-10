library(readr)
library(dplyr)
library(tidyr)
library(kableExtra)
library(DescTools)
# Read akableExtra# Read all results files
# 设置文件夹路径
folder_path <- "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\3T\\step_three"

setwd( "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1")
#### create folder path ####
folder_3T <- data.frame(Model = c("One-step", "Two-step", "Three-step"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\3T\\one_step",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\3T\\step_two",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\3T\\step_three"
                        ))
folder_6T <- data.frame(Model = c("One-step", "Two-step", "Three-step"),
                        Path = c(
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\6T\\one_step",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\6T\\step_two",
                          "C:\\Users\\liuy16\\OneDrive - Universiteit Leiden\\Project_1\\Results\\study1\\1.5\\6T\\step_three"
                        )) # win

folder_3T <- data.frame(Model = c("One-step", "Two-step", "Three-step"),
                        Path = c(
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/3T/one_step",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/3T/step_two",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/3T/step_three"
                        ))
folder_6T <- data.frame(Model = c("One-step", "Two-step", "Three-step"),
                        Path = c(
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/6T/one_step",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/6T/step_two",
                          "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_1/Results/study1/1.5/6T/step_three"
                        )) # mac

folder <- folder_6T


#### read file and complie results of each model ####
library(SimDesign)
sim_conditions <- c("EH1000", "EH500", 
                    "EL1000", "EL500",
                    "UH1000", "UH500",
                    "UL1000", "UL500")
#sim_models <- c("One step A", "One step B", "One step C")

for (i in 1:nrow(folder)) {
  print(i)
  
  folder_path <-folder[i,2]
  # 读取文件夹中所有CSV文件
  results <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE) %>% 
    lapply( read.csv, header = F, sep = ",")
  names(results) <- list.files(folder_path, pattern = "\\.csv$", full.names = F)
  
  model <- folder[i,1]
  
  ##### select data #####
  if (model== "One-step"){
    ###### Create empty matrix 
    ##### set number of simulation and conditions ####
    nsim <- 100
    conditions <- 8
    
    #### Entropy
    Entropy <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results)))
    
    
    ####  parameters  ####
    Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
    Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
    
    #### standard error ####
    SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)

    for (i in 1:length(results)){
      
      # Entropy
      Entropy[ ,i] <- results[[i]][c(seq(3, nrow(results[[ i ]]), by = 10)), 10]
      
      ### Parameters  
      # The intercept of latent class
      Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
      
      # The slope of latent class (covariates X1)
      Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
      
      ### Standard error
      # The intercept of latent class
      SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
      
      # The slope of latent class (covariates X1)
      SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
      
    } 
    
    # Bias
    onestepB_Beta_c <- rbind(
      PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
      AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
      RSE_beta_1c = RSE(SE_1c, Beta_1c)
    ) 
    colnames(onestepB_Beta_c ) <- sim_conditions
    onestep_bias <- onestepB_Beta_c
    
    # CR
    CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
    
    for (j in 1:length(sim_conditions)) {
      CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%")))

                  )
      for(i in 1:nsim) {
        CIs$beta_c1[i, 1] <- Beta_1c[i, j] - 1.96 * SE_1c[i, j]
        CIs$beta_c1[i, 2] <- Beta_1c[i, j] + 1.96 * SE_1c[i, j]
      } # calculate confidence interval for each simulation of conditions
      
      CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
    } # calculate the coverage rates for each condition
    SEP <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
    colnames(SEP) <- sim_conditions 
    onestep_cr <- SEP[-1, ]
    
    # MSE
    MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T)
                 )
    colnames(MSE) <- sim_conditions 
    onestep_mse <- MSE
    
  }
  
 
  if (model== "Two-step"){
    ###### Create empty matrix 
    ##### set number of simulation and conditions ####
    nsim <- 100
    conditions <- 8
    
    ####  parameters  ####
    Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
    Beta_2c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)

    #### standard error ####
    SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)

    for (i in 1:length(results)){
      
      ### Parameters  
      # The intercept of latent class
      Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
      
      # The slope of latent class (covariates X1)
      Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]

      ### Standard error
      # The intercept of latent class
      SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
      
      # The slope of latent class (covariates X1)
      SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
    } 
    # Bias
    onestepB_Beta_c <- rbind(
      PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
      AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
      RSE_beta_1c = RSE(SE_1c, Beta_1c)
    ) 
    colnames(onestepB_Beta_c ) <- sim_conditions
    twostep_bias <- onestepB_Beta_c
    
    # CR
    CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions))
                )
    for (j in 1:length(sim_conditions)) {
      CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%")))
                  )
      for(i in 1:nsim){
        CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
        CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
      } # calculate confidence interval for each simulation of conditions
      CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
  } # calculate the coverage rates for each condition
  # calculate the coverage rates for each condition
    SEP <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
    colnames(SEP) <- sim_conditions 
    twostep_cr <- SEP[-1, ]
    
    # MSE
    MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T)
                 )
    colnames(MSE) <- sim_conditions 
    twostep_mse <- MSE
  }

  
  if (model== "Three-step"){
    ###### Create empty matrix 
    ##### set number of simulation and conditions ####
    nsim <- 100
    conditions <- 8
    
    ####  parameters  ####
    Beta_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    Beta_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)
    
    #### standard error ####
    SE_0c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The intercept of latent class
    SE_1c <- matrix(nrow = nsim, ncol = conditions, dimnames = list(c(1:nsim), names(results))) # The slope of latent class (covariates X1)

    for (i in 1:length(results)){
      
      ### Parameters  
      # The intercept of latent class
      Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 6]
      
      # The slope of latent class (covariates X1)
      Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 7]
      
      ### Standard error
      # The intercept of latent class
      SE_0c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 6]
      
      # The slope of latent class (covariates X1)
      SE_1c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 7]
    }
    # Bias
    onestepB_Beta_c <- rbind(
      PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
      AB_beta_1c = bias(Beta_1c, parameter = -0.5, abs = T),
      RSE_beta_1c = RSE(SE_1c, Beta_1c)
    ) 
    colnames(onestepB_Beta_c ) <- sim_conditions
    threestep_bias <- onestepB_Beta_c
    
    # CR
    CR <- list( beta_c1 = matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions)))
    for (j in 1:length(sim_conditions)) {
      CIs <- list(beta_c1 = matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%")))
                  )
      for(i in 1:nsim){
        CIs$beta_c1[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
        CIs$beta_c1[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
      } # calculate confidence interval for each simulation of conditions
      CR$beta_c1[j] <- ECR(CIs$beta_c1, parameter = c(-0.5))
         } # calculate the coverage rates for each condition
    SEP <- t(BinomCI(CR$beta_c1*100, n=100, conf.level = 0.95, method="wilson"))
    colnames(SEP) <- sim_conditions 
    threestep_cr <- SEP[-1,]
    
    # MSE
    MSE <- rbind(Beta_1c = SimDesign::RMSE(Beta_1c, -0.5, MSE = T)
                 )
    colnames(MSE) <- sim_conditions 
    threestep_mse <- MSE
  }
 }

# beta_1c and beta_2c
#### percentage relative bias ####
table_prb_c1 <- rbind("One-step" = onestep_bias[1, ],
                      "Two-step" = twostep_bias[1, ],
                      "Three-step" = threestep_bias[1, ])

#### absolute bias ####
table_ab_c1 <- rbind("One-step" = onestep_bias[2, ],
                     "Two-step" = twostep_bias[2, ],
                     "Three-step" = threestep_bias[2, ])

#### standard error ####
table_se_c1 <- rbind("One-step" = onestep_bias[3, ],
                     "Two-step" = twostep_bias[3, ],
                     "Three-step" = threestep_bias[3, ])

#### cr ####
table_sep_lw <- rbind("One-step" = onestep_cr[1,],
                     "Two-step" = twostep_cr[1,],
                     "Three-step" = threestep_cr[1,])
table_sep_up <- rbind("One-step" = onestep_cr[2,],
                     "Two-step" = twostep_cr[2,],
                     "Three-step" = threestep_cr[2,])

#### mse ####
table_MSE_c1 <- rbind("One-step" = data.frame(onestep_mse),
                      "Two-step" = data.frame(twostep_mse),
                      "Three-step" = data.frame(threestep_mse)) 



#### detailed table data format ####
results_3t <- list(table_ab_c1=table_ab_c1, table_sep_lw=table_sep_lw,
                   table_sep_up=table_sep_up, table_MSE_c1=table_MSE_c1,
                   table_prb_c1=table_prb_c1, table_se_c1=table_se_c1
                   )

results_6t <- list(table_ab_c1=table_ab_c1, table_sep_lw=table_sep_lw,
                   table_sep_up=table_sep_up, table_MSE_c1=table_MSE_c1,
                   table_prb_c1=table_prb_c1, table_se_c1=table_se_c1
                   )
format_table <- function(Bias, Coverage_rate_lw,Coverage_rate_up) {
  
  
  format_cell <- function(value1, value2, value3) {
    ifelse(value3 >= 0.95,
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
                  caption = "The Absolute Bias ($95\\%$ Confidence Interval of Coverage) Values 
                  Over all Replications for the Latent Class Parameters $\\beta_{x_1}$ and 
                  $\\beta_{x_2}$, in Eight Simulation Conditions \\label{tab:beta_c study1}", 
                  col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
    add_header_above(c(" "=1, "High entropy" = 2, "Low entropy" = 2, "High entropy" = 2, "Low entropy" = 2 )) %>%
    add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
    pack_rows(index=c("$\\beta_{x_1}$ = -0.5" = 3),escape = F) %>% 
    footnote(general = "Two step (Three step) A = Two-step (Three-step) estimator ignore the DE 
             at step 1 model. Two step (Three step) B = Two-step (Three-step) estimator only model 
             the DE on the latent intercept at step 1 model. Two step (Three step) C = Two-step 
             (Three-step) estimator model the DE on both growth factors at step 1 model. " , threeparttable = T) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 7) %>%
    landscape()
  
  print(Table)
  return(formatted_data)
}
mean_results <- results_6t

prb <- mean_results$table_ab_c1
cr_lw <- mean_results$table_sep_lw
cr_up <- mean_results$table_sep_up
se <- mean_results$table_se_c1
bias_stu1 <- format_table(prb, cr_lw, cr_up)
bias_stu16 <- format_table(prb, cr_lw, cr_up)
bias_study1 <- rbind(bias_stu1,bias_stu16)


#### concise table data format ####
Results_3t <- cbind("$Bias$" = rowMeans(table_ab_c1)[1],
                    "$CVR$" = rowMeans(table_cr_c1)[1],
                    "$SE/SD$" = rowMeans(table_se_c1)[1],
                    "$Bias$" = rowMeans(table_ab_c1)[2],
                    "$CVR$" = rowMeans(table_cr_c1)[2],
                    "$SE/SD$" = rowMeans(table_se_c1)[2],
                    "$Bias$" = rowMeans(table_ab_c1)[3],
                    "$CVR$" = rowMeans(table_cr_c1)[3],
                    "$SE/SD$" = rowMeans(table_se_c1)[3])

Results_6t <- cbind("$Bias$" = rowMeans(table_ab_c1)[1],
                    "$CVR$" = rowMeans(table_cr_c1)[1],
                    "$SE/SD$" = rowMeans(table_se_c1)[1],
                    "$Bias$" = rowMeans(table_ab_c1)[2],
                    "$CVR$" = rowMeans(table_cr_c1)[2],
                    "$SE/SD$" = rowMeans(table_se_c1)[2],
                    "$Bias$" = rowMeans(table_ab_c1)[3],
                    "$CVR$" = rowMeans(table_cr_c1)[3],
                    "$SE/SD$" = rowMeans(table_se_c1)[3])

full_results <- rbind(Results_3t, Results_6t); rownames(full_results) <- c("Three time points", "Six time points")

kable(full_results, "latex", digits = 2, booktabs = T, escape =F,
      caption = "The Absolute Bias, Coverage Rates and Standard Error Ratio Values Over All Replications, Averaged Over all Eight Conditions of $\\beta_{1x_1}$ in study 1") %>% 
  add_header_above(c(" "=1,"One-step" = 3, "Two-step" = 3, "Three-step" = 3 ), escape = F) %>%
  #  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("Three time points" = 1, "Six time points" = 1)) %>% 
  footnote(general = "One step = one-step estimator. Two step = two step estimator. Three step = three step estimator." ) %>%
  kable_styling(latex_options = c("hold_position"),font_size = 8) 
#%>%
#  column_spec(1,width = "2cm") 

#mse_signle_results
MSE_1c <- as.matrix(rbind(results_3t$table_MSE_c1, results_6t$table_MSE_c1)) 
rownames(MSE_1c) <- rep(c("One-step", "Two-step", "Three-step"), 2)

kable(MSE_1c, "latex", digits = 3, booktabs = T, escape = F,
      caption = "The Mean Square Error Values Over all Replications for the $\\beta_{x_1]$ in Study 1  \\label{tab:beta_i mse}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  #add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
  pack_rows(index = c("$\\beta_{x_1}=0.50$" = 6 ),escape = F) %>% 
  pack_rows(index = c("Three time points" = 3, "Six time points" = 3)) %>%
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)



# SE_sigle_results
se_1c <- as.matrix(rbind(results_3t$table_se_c1, results_6t$table_se_c1)) # beta_x1
rownames(se_1c) <- rep(c("One-step", "Two-step", "Three-step"), 2)


kable(se_1c, "latex", digits = 2, booktabs = T, escape = F,
      caption = "The Standard Error Ratio Values Over all Replications for the $\\beta_{x_1]$ in Study 1 \\label{tab:beta_i se}", 
      col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 8)), align = "c") %>% 
  add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2,"High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4,"Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("$\\beta_{x_1}=0.50$" = 3 ),escape = F) %>% 
  footnote(general = "One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model the DE on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) model the DE on both growth factors.", threeparttable = T ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 7)

# 读取文件夹中所有CSV文件
results <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE) %>% 
  lapply( read.csv, header = F, sep = ",")
names(results) <- list.files(folder_path, pattern = "\\.csv$", full.names = F)
test_dat <- results$step3_EH_1000_6.csv

###### Create empty matrix 



######## single tables ########
#### set number of simulation and conditions ####
nsim <- 100
conditions <- 8
####  parameters   ####
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

#####

###### One-step 
for (i in 1:length(results)){
  
  # Parameters
  # The intercept of latent class
  Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
  
  # The slope of latent class (covariates X1)
  Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
  
  # The coefficient of latent intercept 
  Beta_0i$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 9]
  Beta_0i$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 10]
  
  # The coefficient for latent slope
  Beta_0s$class1[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 13]
  Beta_0s$class2[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 14]

  
  # standard error
  # The intercept of latent class
  SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
  
  # The slope of latent class (covariates X1)
  SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]
  
} 

###### step two 
for (i in 1:length(results)){
  
  ### Parameters  
  # The intercept of latent class
  Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 6]
  
  # The slope of latent class (covariates X1)
  Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 10)), 7]
  
  ### Standard error
  # The intercept of latent class
  SE_0c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 6]
  
  # The slope of latent class (covariates X1)
  SE_1c[ ,i] <- results[[i]][c(seq(8, nrow(results[[ i ]]), by = 10)), 7]

  
} 

###### Step 3
for (i in 1:length(results)){
  
  ### Parameters  
  # The intercept of latent class
  Beta_0c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 6]
  
  # The slope of latent class (covariates X1)
  Beta_1c[ ,i] <- results[[i]][c(seq(5, nrow(results[[ i ]]), by = 8)), 7]
  
  
  ### Standard error
  # The intercept of latent class
  SE_0c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 6]
  
  # The slope of latent class (covariates X1)
  SE_1c[ ,i] <- results[[i]][c(seq(7, nrow(results[[ i ]]), by = 8)), 7]
  
} 

#####
###### Percentage relative bias(PRB) of parameters #######
library(SimDesign)

sim_conditions <- c("EH1000", "EH500", 
                    "EL1000", "EL500",
                    "UH1000", "UH500",
                    "UL1000", "UL500")
### latent class

onestepB_Beta_c <- rbind(
  PRB_beta_1c = bias(Beta_1c, parameter = -0.5, type = 'relative', percent = T),
  PRB_abs = bias(Beta_1c, parameter = -0.5,abs = T),
  RSE_c1 = RSE(SE_1c, Beta_1c)
 #RSE_beta_2c = RSE(SE_2c, Beta_2c), 
 #PRB_abs = bias(Beta_2c, parameter = -0.5, type = 'relative', percent = T, abs = T)
) 
colnames(onestepB_Beta_c ) <- sim_conditions
onestepB_Beta_c

# storeage results
onestep_bias <- onestepB_Beta_c
steptwo_bias <- onestepB_Beta_c
stepthree_bias <- onestepB_Beta_c

bias_PRC <- rbind("One step" = onestep_bias[1, ],
                    "Two step" = steptwo_bias[1, ],
                    "Three step" = stepthree_bias[1, ])
bias_ABS <- rbind("One step" = onestep_bias[2, ],
                        "Two step" = steptwo_bias[2, ],
                        "Three step" = stepthree_bias[2, ])

table_se <- rbind("One step" = onestep_bias[3, ],
                  "Two step" = steptwo_bias[3, ],
                  "Three step" = stepthree_bias[3, ])
table_se
#### Coverage rates ####
CR <- matrix(NA, 1, 8, dimnames = list(c("CR"), sim_conditions))
for (j in 1:length(sim_conditions)) {
  CIs <- matrix(NA, 100, 2, dimnames = list(rep(1:100), c("2.5%", "97.5%")))
  for(i in 1:nsim){
    CIs[i, 1] <- Beta_1c[i,j] - 1.96*SE_1c[i, j]
    CIs[i, 2] <- Beta_1c[i,j] + 1.96*SE_1c[i, j]
  } # calculate confidence interval for each simulation of conditions
  CR[j] <- ECR(CIs, parameter = c(-0.5))
} # calculate the coverage rates for each condition
CR

# storage results
onestep_cr <- CR
two_step_cr <- CR
three_step_cr <- CR


table_cr <- rbind("One step" = onestep_cr,
                  "Two step" = two_step_cr,
                  "Three step" = three_step_cr)
table_cr

#### MSE ####
MSE <- RMSE(Beta_1c, -0.5, MSE = T)
MSE

# storage results
onestep_mse <- MSE
twostep_mse <- MSE
threestep_mse <- MSE

table_mse <- rbind("One step" = onestep_mse,
                   "Two step" = twostep_mse,
                   "Three step" = threestep_mse)
colnames(table_mse) <- sim_conditions 
table_mse
mse <- rowSums(table_mse)

Study1_results <- rbind("Percentage relative bias" = rowMeans(bias_PRC),
                        "Absolute bias" = rowMeans(bias_ABS),
                        "Mean square error" = rowMeans(table_mse),
                        "Coverage Rates" = rowMeans(table_cr),
                        "Standard error ratio" =rowMeans(table_se)
                        )
Study1_result_3t <- Study1_results
Study1_result_6t <- Study1_results

study1 <- rbind(Study1_result_3t,
                Study1_result_6t)

kable(study1, "latex", digits = 2, booktabs = T, 
      caption = "The Percentahe Relative Bias, Absolute Bias, Mean Square Error, Coverage Rates and Standard Error Ratio Values Over All Replications, Averaged Over all Eight Conditions of Study 1") %>% 
#  add_header_above(c(" "=1, "High entropy" = 2, "Low entropy" = 2, "High entropy" = 2, "Low entropy" = 2 )) %>%
#  add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
  pack_rows(index=c("Three time points" = 5, "Six time points" = 5)) %>% 
  footnote(general = "One step = one-step estimator. Two step = two step estimator. Three step = three step estimator." ) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1,width = "8cm") %>%
  column_spec(4,width = "2cm")




########

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
Long_beta_c_3 <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_3t$table_se_c1) ,var="Models"),
                                          rownames_to_column(as.data.frame(results_3t$table_se_c2) ,var="Models")), 
                                    cols = EH1000:UL500,
                                    names_to = "Conditions",
                                    values_to = "SE_ratio"),
                       Par= "beta_1x1",
                       Time_points = "Three time points")
Long_beta_c_6 <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_6t$table_se_c1) ,var="Models"),
                                          rownames_to_column(as.data.frame(results_6t$table_se_c2) ,var="Models")), 
                                    cols = EH1000:UL500,
                                    names_to = "Conditions",
                                    values_to = "SE_ratio"),
                       Par= "beta_1x1",
                       Time_points = "Six time points")
Long_tab <- rbind(Long_beta_c_3,Long_beta_c_6)


# Define the latex version of factor names
f_names <- list(
  "beta_1x1" = TeX(c("$\\beta_{1x_1}$")),
  "Three time points" = "Three time points" ,
  "Six time points" = "Six time points" 
)
f_labeller <- function (variable, value){
  return(f_names[value])
}

# Generate plots
#gapminder %>% filter(year == 2002) %>%
ggplot(Long_tab, aes(x=SE_ratio,y=Models)) +
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



Long_gamma1_3  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_3t$se_1i_c1), var="Models"),
                                           rownames_to_column(as.data.frame(results_3t$se_1i_c2), var="Models")), 
                                     cols = EH1000:UL500,
                                     names_to = "Conditions",
                                     values_to = "SE_ratio" ),
                        Par = c("gamma_1" ),
                        Time_points = "Three time points")
Long_gamma1_6  <- cbind(pivot_longer(rbind(rownames_to_column(as.data.frame(results_6t$se_1i_c1) ,var="Models"),
                                           rownames_to_column(as.data.frame(results_6t$se_1i_c2) ,var="Models")), 
                                     cols = EH1000:UL500,
                                     names_to = "Conditions",
                                     values_to = "SE_ratio"),
                        Par=c("gamma_1"),
                        Time_points = "Six time points")
Long_tab <- rbind(Long_gamma1_3,Long_gamma1_6)
# Generate plots
#gapminder %>% filter(year == 2002) %>%
ggplot(Long_tab, aes(y=SE_ratio,x=Models)) +
  geom_boxplot() +
  #geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5, fun.args = list(mult = qnorm((1-.95)/2))) +
  #  geom_errorbar(stat="summary", fun.data="mean_se", color="blue", size=1.5) +
  facet_grid(.~Time_points)+
  theme_bw() + coord_flip()+
  labs(y ="Standard Error Ratio", x = "Models")+
  ylim(0.85,1.15)








Entropy_s1_3 <- Entropy
Entropy_s1_6 <-Entropy

Entropy_s1_low <- rbind(Entropy_s1_3[, c(3,4,7,8)], Entropy_s1_6[, c(3,4,7,8)])
Entropy_s1_high <- rbind(Entropy_s1_3[, c(1,2,5,6)], Entropy_s1_6[, c(1,2,5,6)])

summary(Entropy_s1_low)
summary(Entropy_s1_high)






















