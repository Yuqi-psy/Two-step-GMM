
run.template <- function(template.path, envir, temp.filename.base) {
  temp.filename <- paste(temp.filename.base, ".lgs", sep="")
  
  template <- file(template.path, 'r')
  brew(template, output=temp.filename, envir=envir) # writes to file
  close(template)
  
  paste("C:\\Users\\liuy16\\Downloads\\LatentGOLD6.0\\lg60.exe", temp.filename, "/b") #write the correct location ofg LG
} 
