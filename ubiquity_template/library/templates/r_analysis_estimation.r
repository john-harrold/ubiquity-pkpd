#clearing the workspace
rm(list=ls())
# Turning on more verbose error reporting
options(error=traceback)
options(show.error.locations = TRUE)
# Uncomment to set the script directory as the working directory
# This works when calling this file as a script:
# R -e "source('thisfile.r')"
# setwd(dirname(sys.frame(tail(grep('source',sys.calls()),n=1))$ofile))
graphics.off()
library("deSolve")
library("ggplot2")
library("optimx")
library("gdata")
source("library/r_general/ubiquity.r")

# Used for parallelizing 
library("foreach")
library("doParallel")
library("doRNG")

# -------------------------------------------------------------------------

# flowctl       = 'plot previous estimate'
# flowctl       = 'previous estimate as guess'
# flowctl       = 'estimate'
  flowctl       = 'plot guess'
analysis_name   = 'ANAME'
archive_results = TRUE 

# For documentation explaining how to modify the commands below
# See the "R Workflow" section at the link below:
# http://presentation.ubiquity.grok.tv

# -------------------------------------------------------------------------
# Rebuilding the system (R scripts and compiling C code)
build_system(system_file="<SYSTEM_FILE>")

# loading the different functions
source("transient/auto_rcomponents.r")

# Loading the system information
cfg = system_fetch_cfg()
# Initializing the log file ./transient/ubiquity.log
cfg = system_log_init(cfg)
# -------------------------------------------------------------------------

<PSETS>


# # The following will estimate a subset of the parameters:
# pnames = c('PNAME1', 'PNAME2')
# cfg = system_select_set(cfg, "default", pnames)
#
# # This will estiamte all parameters:
# cfg = system_select_set(cfg, "default")

#
# Change initial guess for parameter PNAME to VALUE. The lower bound (lb) and
# upper bound (ub) can also be changed as well
#
# cfg = system_set_guess(cfg, pname="PNAME", value=VALUE, lb=NULL, ub=NULL) 

# -------------------------------------------------------------------------
# Setting options
# 
# Specify output times here using sparse sampling (large time steps) to make
# the estimation quick. See down below where you can specify the sampling to 
# generate smooth profiles when plotting. This will be the
# default output times unless overwritten at the cohort level:
# cfg=system_set_option(cfg, group  = "simulation", 
#                            option = "output_times", 
#                            seq(0,100,1))
#
# cfg=system_set_option(cfg,group = "simulation", option = "solver", value = "lsoda")
#
# To overwrite solver options use the following:
# cfg=system_set_option(cfg,group  = "solver",
#                           option = "atol",   
#                           value  = 1e-10)
#
# cfg=system_set_option(cfg,group  = "solver",
#                           option = "rtol",   
#                           value  = 1e-10)
#
# Uncomment to specify ode file to use. r-file or c-file
# cfg=system_set_option(cfg, group  = "simulation", 
#                            option = "integrate_with",
#                            value  = "r-file")
# 
# Optimizer options
# 
# See the documentation for optim (?optim) for valid values for method and
# elements for the control. You can change the optimizer to optimx and then
# you would just need to adjust the method and control options accordingly.
# The commented code below represents the default values that are used. To
# perform a simulated annealing change the "method" to SANN. 
#
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "optimizer", 
#                              value  = "optim")
#
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "method",
#                              value  = "Nelder-Mead")
# 
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "control", 
#                              value  = list(trace  = TRUE,
#                                            maxit  = 500,
#                                            REPORT = 10))
#
# To use the global optimization routines, uncomment the appropriate lines
# below. IMPORTANT: Because of how the parameter space is sampled it is
# important to set the parameter bounds to reasonable values>
#
# To use particle swarm optimization use the following:
#
# library("pso")
#
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "optimizer", 
#                              value  = "pso")
# 
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "method",
#                              value  = "psoptim")
#
# The following will use the genetic algorithm:
#
# library("GA")
#
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "optimizer", 
#                              value  = "ga")
# 
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "method",
#                              value  = "ga")
# 
# cfg = system_set_option(cfg, group  = "estimation",
#                              option = "control", 
#                              value  = list(maxiter   = 10000, 
#                                            optimArgs = list(method  = "Nelder-Mead",
#                                                             maxiter = 1000)))


# -------------------------------------------------------------------------
# Loading Datasets
#
# From Excel sheet
# cfg = system_load_data(cfg, dsname     = "DSNAME", 
#                             data_file  = "DS.xls", 
#                             data_sheet = "SHEET")
#
# From csv 
# cfg = system_load_data(cfg, dsname     = "DSNAME", 
#                             data_file  = "DS.csv")
#
# From tab 
# cfg = system_load_data(cfg, dsname     = "DSNAME", 
#                             data_file  = "DS.tab")
#    


# -------------------------------------------------------------------------
# Defining the cohorts
#
# Clearing all of the cohorts
cfg = system_clear_cohorts(cfg)
 
# One entry for each cohort:
# For more information type:
#
# help system_define_cohort
#
# It is necessary to replace the following compontents:
#
# CHNAME    - cohort name
# COLNAME   - column name in dataset
# ONAME     - output name
# TIMECOL   - column name in dataset with the observation times
# TS        - model timescale corresponding to TIMECOL
# OBSCOL    - column name in dataset with the observation values
# MODOUTPUT - model output corresponding to OBSCOL
#
# Only specify bolus and infusion inputs that are non-zero. Simply ignore
# those that don't exist for the given cohort. Covariates should be specified
# to overwrite the default covariate values
cohort = c()
cohort$name                                 = 'CHNAME'
cohort$cf$COLNAME                           = c()
cohort$cf$COLNAME                           = c()
cohort$dataset                              = 'DSNAME'

<BOLUS><INFUSION_RATES><COVARIATES>
cohort$outputs$ONAME$of$COLNAME             = c()
cohort$outputs$ONAME$of$COLNAME             = c()

cohort$outputs$ONAME$obs$time               = 'TIMECOL' 
cohort$outputs$ONAME$obs$value              = 'OBSCOL'
cohort$outputs$ONAME$obs$missing            = -1
cohort$outputs$ONAME$model$time             = 'TS'        
cohort$outputs$ONAME$model$value            = 'MODOUTPUT'
cohort$outputs$ONAME$model$variance         = 'PRED^2'
cohort$outputs$ONAME$options$marker_color   = 'black'
cohort$outputs$ONAME$options$marker_shape   = 16
cohort$outputs$ONAME$options$marker_line    = 1 

cfg = system_define_cohort(cfg, cohort)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Defining the cohorts with dosing from a NONMEM file
# 
#  DS - name of dataset defined using system_load_dataset
# 
# Include all records in the dataset
#  filter = NULL
# 
# Include only records matching the following filter
#  filter = list()
#  filter$COLNAME = c()
# 
# Mapping information: 
# 
# Inputs:
#  INPUTMAP = list()
#  INPUTMAP$bolus$SPECIES$CMT_NUM            =  1
#  INPUTMAP$infusion_rates$RATE$CMT_NUM      =  1
#  INPUTMAP$covariates$CVNAME$col_COV        = 'CNAME'
#
# Outputs
#  OBSMAP = list()
#  OBSMAP$ONAME=list(variance     = 'PRED^2',
#                    CMT          =  1,
#                    output       = '<O>',
#                    missing      =  NULL )
#
#  cfg = system_define_cohorts_nm(cfg        =  cfg, 
#                                  filter    =  filter,
#                                  INPUTS    =  INPUTMAP,
#                                  OBS       =  OBSMAP,
#                                  col_GROUP =  NULL,  
#                                  group     =  FALSE)
#
# To overwrite column name mapping the following can be used
#  col_ID    - unique subject identifier 
#  col_CMT   - compartment
#  col_DV    - observations
#  col_TIME  - system time
#  col_AMT   - infusion/dose amounts
#               - these need to be in the same units specified in the system.txt file
#  col_RATE  - rate of infusion or . for bolus
#  col_EVID  - evid (0 - observation, 1 dose)
#  col_GROUP - used to group subjects for plotting purposes. For example if
#              the cohorts were grouped by dose and the column DOSE contained
#              this information, then this will be set to 'DOSE'.
# 
#  cfg = system_define_cohorts_nm( cfg       = cfg, 
#                                  DS        = 'DSNAME',
#                                  col_ID    = 'ID',
#                                  col_CMT   = 'CMT',
#                                  col_DV    = 'DV',
#                                  col_TIME  = 'TIME',
#                                  col_AMT   = 'AMT',
#                                  col_RATE  = 'RATE',
#                                  col_EVID  = 'EVID',
#                                  col_GROUP =  NULL,  
#                                  filter    =  filter,
#                                  INPUTS    =  INPUTMAP,
#                                  OBS       =  OBSMAP,
#                                  group     =  FALSE)

# -------------------------------------------------------------------------
# performing estimation or loading guess/previous results
pest = system_estimate_parameters(cfg, 
                                  flowctl         = flowctl, 
                                  analysis_name   = analysis_name, 
                                  archive_results = archive_results)

# -------------------------------------------------------------------------

# Here you can specify the output times used for the simulations that will
# be used for plotting. Here you want frequent sampling (small time steps)
# so that you have smoother profiles.
# cfg=system_set_option(cfg, group  = "simulation", 
#                            option = "output_times", 
#                            seq(0,100,1))

# Simulating the system at the estimates
erp = system_simulate_estimation_results(pest = pest, cfg = cfg) 
# -------------------------------------------------------------------------


plot_opts = c()
# To customize the figures
# Replace ONAME with the specific outputs defined in the cohorts above. Then
# change the values
# plot_opts$outputs$ONAME$yscale   = 'linear' # 'linear' or 'log'
# plot_opts$outputs$ONAME$ylim     = c(0,1)   # NULL
# plot_opts$outputs$ONAME$xlim     = c(0,1)   # NULL
# plot_opts$outputs$ONAME$ylabel   = ''       # output ONAME
# plot_opts$outputs$ONAME$xlabel   = ''       # cohort time units
# plot_opts$tc$width               = 10       # width  of timecourse plot
# plot_opts$tc$height              = 5.5      # height of timecourse plot
# plot_opts$op$width               = 10       # width  of observed vs predicted
# plot_opts$op$height              = 8.0      # height of observed vs predicted

# Plotting the simulated results at the estimates 
# These figures will be placed in output/
system_plot_cohorts(erp, plot_opts, cfg, prefix=analysis_name)
