#clearing the workspace
rm(list=ls())
# Turning on more verbose error reporting
options(error=traceback)
options(show.error.locations = TRUE)


required_packages= c("deSolve",
                     "ggplot2",
                     "foreach",
                     "shiny",
                     "rhandsontable",
                     "rstudioapi",
                     "optimx",
                     "doParallel",
                     "doRNG")



for(pkg in required_packages){

   # If the package insn't installed we install it
   if(!is.element(pkg, installed.packages()[,1])){
     cat(sprintf(" Package not found, attempting to install: %s ", pkg))
     install.packages(pkg)
   }

}
