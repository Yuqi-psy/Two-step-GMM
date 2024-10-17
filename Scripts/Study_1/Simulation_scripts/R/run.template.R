
run.template <- function(template.path, envir, temp.filename.base) {
  temp.filename <- paste(temp.filename.base, ".lgs", sep="")
  
  template <- file(template.path, 'r')
  brew(template, output=temp.filename, envir=envir) # writes to file
  close(template)
  
  paste("~/lg60.exe", temp.filename, "/b") # specify the location of LatentGOLD
} 
