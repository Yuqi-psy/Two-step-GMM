library(knitr)
library(kableExtra)

# Function to create LaTeX tables with \tbl{} and \tabnote{}
create_latex_table <- function(data, type, par, study, footnote = NULL, label = NULL) {
  
  if(study == "Study 1"){
    
  }
  
  if (type == "Bias"){
    kable_table <-  kable(data, "latex", booktabs = T, escape = F,
                          col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
      add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2 )) %>%
      add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
      pack_rows(index=c("$\\beta_{x_1}$ = -0.5" = 9, "$\\beta_{x_2}$ = 0.75" = 9),escape = F) %>% 
      kable_styling(latex_options = c("hold_position"), position = "center")
    
    # Split the table output to isolate the tabular part
    table_lines <- strsplit(kable_table, "\n")[[1]]
    
    # Escape backslashes in the patterns
    tabular_lines <- table_lines[grep("\\\\begin\\{tabular\\}", table_lines):grep("\\\\end\\{tabular\\}", table_lines)]
    
    # Combine the isolated lines back into a single string
    tabular_content <- paste(tabular_lines, collapse = "\n")
    
    
    # Wrap with \tbl{} and add footnote if provided
    final_table <- paste0(
      "\\tbl{", caption, "}\n",
      "{",   tabular_content, "}"
    )
  }
  
  if (type == "MSE"){
    
    criterias <- "Mean Square Error"
    decimal <- 3
  }
  
  if (type == "SE"){
    
    criterias <- "Standard Error Ratio"
    decimal <- 3
  }
  
  
    if (par == "beta_x1x2"){
      kable_table <- kable(data, "latex", digits = decimal, booktabs = T, escape = F,
                           col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
        add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
        add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
        #add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
        pack_rows(index = c("Three time points" = 18, "Six time points" = 18)) %>%
        pack_rows(index=c("$\\beta_{x_1}=0.50$" = 9, "$\\beta_{x_2}=0.75$" = 9, "$\\beta_{x_1}=0.50$" = 9, "$\\beta_{x_2}=0.75$" = 9),escape =F) %>% 
        kable_styling(latex_options = c("hold_position"), position = "center")
      
      caption <- paste0("The ", criterias, "Values of the Class Membership Coefficients $\\beta_{x_1}$ and $\\beta_{x_2}$ Over 100 Replications for ", study, ".")
      
      footnote <- "\\textit{Note. } One-step A, One-step B, and One-step C are the one-step estimators
                              ignoring the direct effects (specifications A), specifying the direct effect on the latent intercept (specification B), 
                              and on both growth factors (specification C), respectively. Two-step A, Two-step B, and Two-step C are the proposed 
                              two-step estimators with specifications A, B, and C, respectively. Three-step A, Three-step B, and Three-step C are 
                              the three-step estimators with specifications A, B, and C, respectively. $N$ is the total sample size. The bold numbers 
                              reflect conditions that do not contain the 95\\% confidence interval."
      label <- paste0("tab:", study, par, type)
    
    }
    
    if (par == "gamma0"){
      kable_table <- kable(data, "latex", digits = decimal, booktabs = T, escape = F,
                           col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
        add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
        add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
        #add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
        pack_rows(index = c("Three time points" = 8, "Six time points" = 8)) %>%
        pack_rows(index=c("$\\gamma_0^{1}=0.50$" = 4, "$\\gamma_0^{2}=-0.50$" = 4, "$\\gamma_0^{1}=0.50$" = 4, "$\\gamma_0^{2}=-0.50$" = 4),escape =F) %>% 
        kable_styling(latex_options = c("hold_position"), position = "center")
      
      caption <- paste0("The ", criterias, " Values of the Latent Intercept Coefficients $\\gamma_0$ Over 100 Replications for ", study, ".")
      
      footnote <- "\\textit{Note.}One-step B is the one-step estimator that only models the direct effect (DE) 
                              on the latent intercept (specification B). One-step C is the one-step estimator that models 
                              the DE on both growth factors (specification C). Step-1 B is the step-1 model of the proposed 
                              two-step estimator and the three-step estimator with specification B.  Step-1 C is the step-1 
                              model of the proposed two-step estimator and the three-step estimator with specification C. 
                              $\\gamma_0^{(1)}$ is the regression coefficient of latent intercept in Class 1. $\\gamma_0^{(2)}$ 
                                is the regression coefficient of latent intercept in Class 2. $N$ is the total sample size. 
                              The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval."
      label <- paste0("tab:", study, par, type)
    }
    
    if (par == "gamma1"){
      kable_table <- kable(data, "latex", digits = decimal, booktabs = T, escape = F, 
                           col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
        add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
        add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
        #add_header_above(c(" "=1, "Three time points" = 4, "Six time points" = 4)) %>%
        pack_rows(index = c("Three time points" = 4, "Six time points" = 4)) %>%
        pack_rows(index=c("$\\gamma_1^{1}=0.25$" = 2, "$\\gamma_1^{2}=-0.25$" = 2, "$\\gamma_1^{1}=0.25$" = 2, "$\\gamma_1^{2}=-0.25$" = 2),escape =F) %>% 
        kable_styling(latex_options = c("hold_position"), position = "center")
       
      caption <- paste0("The ", criterias, " Values of the Latent Slope Coefficients $\\gamma_1$ Over 100 Replications for ", study, ".")
      
      footnote <- "\\textit{Note.} One-step C is the one-step estimator that models 
                              the DE on both growth factors (specification C).Step-1 C is the step-1 
                              model of the proposed two-step estimator and the three-step estimator with specification C. 
                              $\\gamma_1^{(1)}$ is the regression coefficient of latent slope in Class 1. $\\gamma_1^{(2)}$ 
                                is the regression coefficient of latent slope in Class 2. $N$ is the total sample size. 
                              The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval."
      label <- paste0("tab:", study, par, type)
    }
  
  # Split the table output to isolate the tabular part
  table_lines <- strsplit(kable_table, "\n")[[1]]
  
  # Escape backslashes in the patterns
  tabular_lines <- table_lines[grep("\\\\begin\\{tabular\\}", table_lines):grep("\\\\end\\{tabular\\}", table_lines)]
  
  # Combine the isolated lines back into a single string
  tabular_content <- paste(tabular_lines, collapse = "\n")
  
  
  # Wrap with \tbl{} and add footnote if provided
  final_table <- paste0(
    "\\tbl{", caption, "}\n",
    "{",   tabular_content, "}"
  )
  
  if (!is.null(footnote)) {
    final_table <- paste0(final_table, "\n\\tabnote{", footnote, "}")
  }
  
  if (!is.null(label)) {
    final_table <- paste0("\\begin{table}\n", final_table, "\n\\label{", label, "}\n\\end{table}\n")
  }
  
  return(final_table)
}

create_latex_table_stu1 <- function(data, type, footnote = NULL, label = NULL) {
  
  if (type == "BIAS"){
    
    criterias <- "Absolute Bias ($95\\%$ Confidence Interval of Coverage Rates)"
    decimal <- 3
  }
  
  if (type == "MSE"){
    
    criterias <- "Mean Square Error"
    decimal <- 3
  }
  
  if (type == "SE"){
    
    criterias <- "Standard Error Ratio"
    decimal <- 2
  }

  
  kable_table <- kable(data, "latex", digits = decimal, booktabs = T, escape = F,
         col.names = c(rep(c("$N$ = 1000", "$N$ = 500"), 4)), align = "c") %>% 
      add_header_above(c(" "=1, "High entropy" = 2, "Moderate entropy" = 2, "High entropy" = 2, "Moderate entropy" = 2)) %>%
      add_header_above(c(" "=1, "Mixing ratio = 0.50/0.50" = 4, "Mixing ratio = 0.30/0.70" = 4)) %>%
      #add_header_above(c(" "=1, "Three time points" = 8, "Six time points" = 8)) %>%
      pack_rows(index = c("$\\beta_{x_1}=0.50$" = 6 ),escape = F) %>% 
      pack_rows(index = c("Three time points" = 3, "Six time points" = 3)) %>%
      kable_styling(latex_options = c("hold_position"), position = "center")

    
    caption <- paste0("The ", criterias, " Values of the Class Membership Coefficient $\\beta_{x_1}$ Over 100 Replications for Study 1.")
    
    footnote <- "\\textit{Note. } One step (Step one) B = one-step estimator (step one model of step-wise estimator) only model 
               the direct effect (DE) on the latent intercept. One step (step one) C = one-step estimator (step one model of step-wise estimator) 
               model the DE on both growth factors."
    label <- paste0("tab: study 1 beta_x1", type)
    
  # Split the table output to isolate the tabular part
  table_lines <- strsplit(kable_table, "\n")[[1]]
  
  # Escape backslashes in the patterns
  tabular_lines <- table_lines[grep("\\\\begin\\{tabular\\}", table_lines):grep("\\\\end\\{tabular\\}", table_lines)]
  
  # Combine the isolated lines back into a single string
  tabular_content <- paste(tabular_lines, collapse = "\n")
  
  
  # Wrap with \tbl{} and add footnote if provided
  final_table <- paste0(
    "\\tbl{", caption, "}\n",
    "{",   tabular_content, "}"
  )
  
  if (!is.null(footnote)) {
    final_table <- paste0(final_table, "\n\\tabnote{", footnote, "}")
  }
  
  if (!is.null(label)) {
    final_table <- paste0("\\begin{table}\n", final_table, "\n\\label{", label, "}\n\\end{table}\n")
  }
  
  return(final_table)
}

# Generate LaTeX table
latex_table <- create_latex_table(data = MSE_1i,
                                  type = "MSE",
                                  par = "gamma0",
                                  study = "Study 3",
                                  footnote = T,
                                  label = T)
cat( latex_table)
writeClipboard(latex_table)

latex_table <- create_latex_table_stu1(data = bias_study1,
                                  type = "BIAS",
                                  footnote = T,
                                  label = T)
cat( latex_table)
writeClipboard(latex_table)


latex_table <- create_latex_table(data = MSE_1s,
                                  type = "MSE",
                                  par = "gamma1",
                                  caption =  "The Mean Square Error Values Over all Replications for the $\\gamma_1$ in Eight Simulation Conditions for study 2.", 
                                  footnote = "\\textit{Note.}One step (Step one) B = one-step estimator (step one model of step-wise estimator) 
           only model the DE on the latent intercept. One step (step one) C = one-step estimator 
           (step one model of step-wise estimator) model the DE on both growth factors.")
cat("\\begin{table}\n", latex_table, "\n\\label{tab:gamma_1 mse}\n\\end{table}\n")

latex_table <- create_latex_table(data = se_1c_long,
                                  type = "SE",
                                  par = "beta_x1x2",
                                  caption =  "The Standard Error Ratio Values Over all Replications for the Latent Class Parameters $\\gamma_0$ in Study 2.", 
                                  footnote = "\\textit{Note.}One-step B is the one-step estimator that only models the direct effect (DE) 
             on the latent intercept (specification B). One-step C is the one-step estimator that models 
             the DE on both growth factors (specification C). Step-1 B is the step-1 model of the proposed 
             two-step estimator and the three-step estimator with specification B.  Step-1 C is the step-1 
             model of the proposed two-step estimator and the three-step estimator with specification C. 
             $\\gamma_0^{(1)}$ is the regression coefficients of latent intercept in Class 1. $\\gamma_0^{(2)}$ 
             is the regression coefficients of latent intercept in Class 2. $N$ is the total sample size. 
             The bold numbers reflect the nominal level of coverage rates not within the 95\\% confidence interval.")
cat("\\begin{table}\n", latex_table, "\n\\label{tab:beta_i se_sig_stu3}\n\\end{table}\n")