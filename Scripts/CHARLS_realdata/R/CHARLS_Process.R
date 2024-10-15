# Filtering individuals with diabetes
blood_CHARLS <- read_dta(" ") # loading baseline data with individual blood information
filtered_blood <- blood_CHARLS[blood_CHARLS$newhba1c >= 6.5 |
                                 blood_CHARLS$newglu >= 126, ] # selecting individuals according to their fasting plasma glucose & glycosylated hemoglobin
ID_blood <- filtered_blood[!is.na(filtered_blood$ID), ]
filtered_ids <- ID_blood$ID # extracting baseline diabetes ID of selected individuals
raw_CHARLS <- read_dta(" ") # loading four waves data
raw_CHARLS <- raw_CHARLS[!is.na(raw_CHARLS$ID), ]
CHARLS_diab <- raw_CHARLS[raw_CHARLS$ID_w1 %in% filtered_ids, ]
cohort_CHARLS <- CHARLS_diab[which(CHARLS_diab$hacohort_c == 1), ]
summary(cohort_CHARLS$hacohort_c) # checking the sample cohort 

# Selecting varaibels
CHARLS_cog <- cohort_CHARLS %>% select( 
  
  # Measures-Executive function
  r1ser7, r2ser7, r3ser7, r4ser7,#  serial 7's
  r1orient,r2orient, r3orient, r4orient, # summary of date naming
  r1draw, r2draw, r3draw, r4draw, # drawing picture
  
  # Measures-Episodic memory
  r1tr20, r2tr20, r3tr20, r4tr20,# recall summary score(immdiate recall &delay recall)
  
  # covariates time invariant
  rabyear, # age
  ragender, # gender
  "raeducl"= raeduc_c, # educational attachment
  h1rural, h2rural, h3rural, h4rural, # level of urbanlization
  r1cesd10, r2cesd10, r3cesd10, r4cesd10,
  r1cesd10m, r2cesd10m, r3cesd10m, r4cesd10m # depression score
) 


# Re-code age variable
categorial_age <- c(-Inf, 1949, 1964, Inf)
labels <- c("very old", "old", "young")
CHARLS_cog$age <- cut(CHARLS_cog$rabyear, categorial_age, labels, ordered_result = T)
CHARLS_cog <- CHARLS_cog[!is.na(CHARLS_cog$age), ]

# Re-code education variable
categorial_edu <- c(-Inf, 1, 4, Inf)
labs <- c("no formal education", "primary school", "middle school and above")
CHARLS_cog$raeducl <- cut(CHARLS_cog$raeducl, categorial_edu, labs, ordered_result = T)
CHARLS_cog$ragender <- as.character(CHARLS_cog$ragender)

# Re-code gender variable
CHARLS_cog <- CHARLS_cog %>%
  mutate(ragender = recode(ragender, "1" = "male", "2" = "female"))

# Compute average of recalling score
CHARLS_cog$r1tr20 <- CHARLS_cog$r1tr20 / 2
CHARLS_cog$r2tr20 <- CHARLS_cog$r2tr20 / 2
CHARLS_cog$r3tr20 <- CHARLS_cog$r3tr20 / 2
CHARLS_cog$r4tr20 <- CHARLS_cog$r4tr20 / 2

# Compute sum of score of mental intactness
CHARLS_cog$mental_intact1 <- CHARLS_cog$r1ser7 + CHARLS_cog$r1orient + CHARLS_cog$r1draw
CHARLS_cog$mental_intact2 <- CHARLS_cog$r2ser7 + CHARLS_cog$r2orient + CHARLS_cog$r2draw
CHARLS_cog$mental_intact3 <- CHARLS_cog$r3ser7 + CHARLS_cog$r3orient + CHARLS_cog$r3draw
CHARLS_cog$mental_intact4 <- CHARLS_cog$r4ser7 + CHARLS_cog$r4orient + CHARLS_cog$r4draw

# Compute sum scores of cognitive function
CHARLS_cog$cog_function1 <- CHARLS_cog$mental_intact1 + CHARLS_cog$r1tr20
CHARLS_cog$cog_function2 <- CHARLS_cog$mental_intact2 + CHARLS_cog$r2tr20
CHARLS_cog$cog_function3 <- CHARLS_cog$mental_intact3 + CHARLS_cog$r3tr20
CHARLS_cog$cog_function4 <- CHARLS_cog$mental_intact4 + CHARLS_cog$r4tr20

# Exclude participants with repeated measures less than 2 points
data_filterd <- CHARLS_cog %>% rowwise() %>%
  mutate(non_missing = sum(!is.na(c_across(
    cog_function1:cog_function4
  )))) %>%
  filter(non_missing >= 2) %>%
  select(-non_missing)


# Add id variable
id <- c(1:nrow(data_filterd))
CHARLS_cog <- cbind(data_filterd, id)

# Write csv file
write.table(
  CHARLS_cog,
  "~/CHARLS_cog.csv",
  row.names = FALSE,
  col.names = TRUE,
  sep = ",",
  quote = FALSE
)
