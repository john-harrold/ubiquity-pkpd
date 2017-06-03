#clearing the workspace
rm(list=ls())
# Turning on more verbose error reporting
options(error=traceback)
options(show.error.locations = TRUE)
# Uncomment to set the script directory as the working directory
# setwd(dirname(sys.frame(tail(grep('source',sys.calls()),n=1))$ofile))
graphics.off()
library("deSolve")
library("ggplot2")
require("gdata")
source("library/r_general/ubiquity.r");

# Used for parallelizing 
library("foreach")
library("doParallel")
library("doRNG")

# For documentation explaining how to modify the commands below
# See the "R Workflow" section at the link below:
# http://presentation.ubiquity.grok.tv

# Rebuilding the system (R scripts and compiling C code)
build_system(system_file="<SYSTEM_FILE>")

# loading the different functions
source("transient/auto_rcomponents.r");

# Loading the system information
cfg = system_fetch_cfg()
<PSETS>
cfg = system_select_set(cfg, "default")

# fetching the parameter values
parameters = system_fetch_parameters(cfg)

# The previous statement sets 'parameters' to the values 
# in the currently selected parameter set. To overwrite 
# a specific parameter uncomment the following statement 
# and replace PNAME with the name of the parameter 
# and VALUE with the desired value:
#
# parameters$PNAME = VALUE;

<BOLUS>
<INFUSION_RATES>
<COVARIATES>
<OUTPUT_TIMES>


# The following applies to both individual and stochastic simulations:
# Define the solver to use
cfg=system_set_option(cfg,group = "simulation", option = "solver", value = "lsoda")
#
# To overwrite solver options use the following:
# cfg=system_set_option(cfg,group  = "solver",
#                           option = "atol",   
#                           value  = 1e-10)
# cfg=system_set_option(cfg,group  = "solver",
#                           option = "rtol",   
#                           value  = 1e-10)

# Specify the output times 
cfg=system_set_option(cfg, group  = "simulation", 
                           option = "output_times", 
                           seq(0,100,1))
# By default, important times will be included in the simulation output
# e.g., bolus times, sampling before and after rate switches, etc.
# uncomment to only evalutate at the output times specified above
# cfg=system_set_option(cfg, group  = "simulation", 
#                            option = "include_important_output_times", 
#                            value  = "no")
# Uncomment to specify ode file to use
# cfg=system_set_option(cfg, group  = "simulation", 
#                            option = "integrate_with",
#                            value  = "r-file")

# -------------------------------------------------------------------------
# Individual Simulation:
som = run_simulation_ubiquity(parameters, cfg)
# # replace TS     with a timescale (i.e. days) and 
# #         OUTPUT with a named output  (i.e. Cp)
#plot(som$simout$TS,        som$simout$OUTPUT)
# myfig = ggplot() + 
#         geom_line(data=som$simout, aes(x=ts.TS,   y=OUTPUT), color="red") 
# print(myfig)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# # Stochastic Simulation:
# # To use this you need to have specified variability 
# # in your model using the <IIV:?>, <IIV:?:?>, and <IIVCOR:?:?> 
# # delimiters see examples/system-iiv.txt
#
# cfg = system_set_option(cfg, group="stochastic", option="nsub",    value=100)
# cfg = system_set_option(cfg, group="stochastic", option="ci",      value=95 )
# cfg = system_set_option(cfg, group="stochastic", option="seed",    value=8675309)
# cfg = system_set_option(cfg, group="stochastic", option="ponly",   value=FALSE)
# cfg = system_set_option(cfg, group="stochastic", option="states",  value=list())
# cfg = system_set_option(cfg, group="stochastic", option="outputs", value=c('OP1', 'OP2'))
#
# som   = simulate_subjects(parameters, cfg)
# 
# # To parallelize the simulations uncomment the following:
#  cfg=system_set_option(cfg, group  = "simulation",
#                             option = "parallel",    
#                             value  = "multicore")
#  
#  cfg=system_set_option(cfg, group  = "simulation",
#                             option = "compute_cores", 
#                             value  = detectCores() - 1)
#  
# #
# # replace TS     with a timescale (i.e. days) and 
# #         OUTPUT with a named output  (i.e. Cp)
#
# p = ggplot() 
# p = p + geom_ribbon(data=som$tcsummary, aes(x=ts.TS, ymin=o.OUTPUT.lb_ci, ymax=o.OUTPUT.ub_ci), fill='cadetblue1', alpha=0.6)
# p = p + geom_line(  data=som$tcsummary, aes(x=ts.TS, y=o.OUTPUT.median, color='OUTPUT'), linetype='solid', size=0.9) 
# p = p + xlab('Time (TS)')                     
# p = p + ylab('Output')                        
# p = prepare_figure('present', p)              
# p = p + scale_colour_manual(values=c("OUTPUT"="darkblue"))   
# p = p + theme(legend.title = element_blank()) 
# p = p + theme(legend.position = 'bottom') 
#
# print(p)
# 
# -------------------------------------------------------------------------
