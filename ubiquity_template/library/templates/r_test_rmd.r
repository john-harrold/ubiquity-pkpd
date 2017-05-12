

Rmdfile = "my_report.Rmd"
load("transient/rgui/default/gui_som.RData")
load("transient/rgui/default/gui_state.RData")
params = list()
params$cfg = cfg
params$som = som
rmarkdown::render(Rmdfile,
                  params        = params,
                  output_format = "html_document")
