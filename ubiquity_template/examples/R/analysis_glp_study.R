#clearing the workspace
rm(list=ls())
graphics.off()
options(show.error.locations = TRUE)

# If we cannot load the ubiquity package we try the stand alone distribution
if("ubiquity" %in% rownames(installed.packages())){require(ubiquity)} else 
{source(file.path('library', 'r_general', 'ubiquity.R')) }

# For documentation explaining how to modify the commands below
# See the "R Workflow" section at the link below:
# http://presentation.ubiquity.grok.tv

# Rebuilding the system (R scripts and compiling C code)
cfg = build_system(system_file="system.txt",
                   output_directory     = file.path(".", "output"),
                   temporary_directory  = file.path(".", "transient"))


# set name                  | Description
# -------------------------------------------------------
# default                   | NHP
# human                     | Human


# fetching the nhp and human parmeters:
cfg         = system_select_set(cfg, "default")
p_nhp       = system_fetch_parameters(cfg)
cfg         = system_select_set(cfg, "human")
p_human     = system_fetch_parameters(cfg)

# Creating an empty report
cfg = system_report_init(cfg)
cfg = system_report_slide_title(cfg,
        title   = "GLP Tox Study Design")







# Saving the report
system_report_save(cfg, output_file=file.path("output", "glp_study_report.pptx"))
