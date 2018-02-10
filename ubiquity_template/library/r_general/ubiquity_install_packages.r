#clearing the workspace
rm(list=ls())
# Turning on more verbose error reporting
options(error=traceback)
options(show.error.locations = TRUE)


required_packages= c("deSolve",
                     "ggplot2",
                     "foreach",
                     "shiny",
                     "gdata",
                     "rmarkdown",     
                     "rhandsontable",
                     "rstudioapi",
                     "optimx",
                     "doParallel",
                     "doRNG")


reccomended_packages = c( "GGally",
                          "officer",
                          "flextable",
                          "gridExtra",
                          "gridGraphics")

for(pkg in required_packages){

   # If the package insn't installed we install it
   if(!is.element(pkg, installed.packages()[,1])){
     cat(sprintf(" Package not found, attempting to install: %s ", pkg))
     install.packages(pkg)
   }

}


# pkgs = search()

# install.packages("devtools")
# library("devtools")
# require("devtools")
# devtools::install_github("AnalytixWare/ShinySky")
