
# Set working directory
setwd(" ")

library(brew)
source(
  "C:~\\run.template.R"
)

# Class enumeration process
envir <- new.env()
brewfile1 <- paste("Classenum_CHARLS", ".brew", sep = "")
run.template(brewfile1, envir = envir, temp.filename.base = "Classenum_CHARLS")
lgsfile1 <- paste("Classenum_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile1,
  "/b"
))	

# Specifying the optimal number of classes
classnumber  <- 2

# Label the specification types
modeltype <- c("A", "B", "C")

# One-step full information ML model
envir <- new.env() 
assign("classnumber", classnumber, envir = envir) # label the optimal number of classes
brewfile2 <- paste("OnestepA_CHARLS", ".brew", sep = "")
run.template(brewfile2, envir = envir, temp.filename.base = "OnestepA_CHARLS")
lgsfile2 <- paste("OnestepA_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile2,
  "/b"
)) # Specification A

envir <- new.env() 
assign("classnumber", classnumber, envir = envir)
brewfile3 <- paste("OnestepB_CHARLS", ".brew", sep = "")
run.template(brewfile3, envir = envir, temp.filename.base = "OnestepB_CHARLS")
lgsfile3 <- paste("OnestepB_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile3,
  "/b"
)) # Specification B

envir <- new.env()
assign("classnumber", classnumber, envir = envir)
brewfile4 <- paste("OnestepC_CHARLS", ".brew", sep = "")
run.template(brewfil4, envir = envir, temp.filename.base = "OnestepC_CHARLS")
lgsfile4 <- paste("OnestepC_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile4,
  "/b"
))	# Specification C 


# Step-wise estimators
envir <- new.env() 
assign("classnumber", classnumber, envir = envir)
assign("modeltype", modeltype[1], envir = envir) # label the specification of direct effects
brewfile5 <- paste("Step1A_CHARLS", ".brew", sep = "")
run.template(brewfile5, envir = envir, temp.filename.base = "Step1A_CHARLS")
lgsfile5 <- paste("Step1A_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile5,
  "/b"
)) # step-one model with specification A


DAT <- paste("step2pars_CHARLS_", "A", ".dat", sep = "")
bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
envir <- new.env()
assign("classnumber", classnumber, envir = envir)

assign("step1a_pars", bpars[-1, ], envir = envir)
brewfile6 <- paste("Step2A_CHARLS", ".brew", sep = "")
run.template(brewfile6, envir = envir, temp.filename.base = "Step2A_CHARLS")
lgsfile6 <- paste("Step2A_CHARLS" , ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile6,
  "/b"
)) # step-two model	 

envir <- new.env() 
assign("classnumber", classnumber, envir = envir)
assign("modeltype", modeltype[2], envir = envir)
brewfile7 <- paste("Step1B_CHARLS", ".brew", sep = "")
run.template(brewfile7, envir = envir, temp.filename.base = "Step1B_CHARLS")
lgsfile7 <- paste("Step1B_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile7,
  "/b"
)) 	# step-one model with specification B 

DAT <- paste("step2pars_CHARLS_", "B", ".dat", sep = "")
bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
envir <- new.env()
assign("classnumber", classnumber, envir = envir)

assign("step1a_pars", bpars[-c(1:3), ], envir = envir)
brewfile8 <- paste("Step2B_CHARLS", ".brew", sep = "")
run.template(brewfile8, envir = envir, temp.filename.base = "Step2B_CHARLS")
lgsfile8 <- paste("Step2B_CHARLS" , ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile8,
  "/b"
)) # step-two model	 

envir <- new.env() 
assign("classnumber", classnumber, envir = envir)
assign("modeltype", modeltype[3], envir = envir)
brewfile9 <- paste("Step1C_CHARLS", ".brew", sep = "")
run.template(brewfile9, envir = envir, temp.filename.base = "Step1C_CHARLS")
lgsfile9 <- paste("Step1C_CHARLS", ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile9,
  "/b"
)) # step-one model with specification C	

DAT <- paste("step2pars_CHARLS_", "C", ".dat", sep = "")
bpars <- as.matrix(read.table(DAT, sep = " ", dec = ","))
envir <- new.env()
assign("classnumber", classnumber, envir = envir)

assign("step1a_pars", bpars[-c(1:3), ], envir = envir)
brewfile10 <- paste("Step2C_CHARLS", ".brew", sep = "")
run.template(brewfile10, envir = envir, temp.filename.base = "Step2C_CHARLS")
lgsfile10 <- paste("Step2C_CHARLS" , ".lgs", sep = "")
shell(paste(
  "C:~\\lg60.exe", # the latentGOLD address
  lgsfile10,
  "/b"
))	# step-two model


envir <- new.env()
assign("modeltype", modeltype[1], envir = envir)
brewfile11 <- paste("Step3A_CHARLS", ".brew", sep = "")
run.template(brewfile11,
             envir = envir,
             temp.filename.base = "Step3_CHARLS")
lgsfile11 <- paste("Step3_CHARLS", ".lgs", sep = "")
shell(
  paste(
    "C:~\\lg60.exe", # the latentGOLD address
    lgsfile11,
    "/b"
  )
) # step-three model (specification A)


for (i in modeltype[-1]) {
  # step 3 model
  envir <- new.env()
  assign("modeltype", i, envir = envir)
  assign("model", i, envir = envir)
  brewfile12 <- paste("Step3_CHARLS", ".brew", sep = "")
  run.template(brewfile12,
               envir = envir,
               temp.filename.base = "Step3_CHARLS")
  lgsfile12 <- paste("Step3_CHARLS", ".lgs", sep = "")
  shell(
    paste(
      "C:~\\lg60.exe", # the latentGOLD address
      lgsfile12,
      "/b"
    )
  )	
}
