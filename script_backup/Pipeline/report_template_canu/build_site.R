#! /usr/bin/env R
#Set our working directory. 
#This helps avoid confusion if our working directory is 
#not our site because of other projects we were 
#working on at the time. 

library(rmarkdown)
setwd("report_template_canu")

#render your sweet site. 
#rmarkdown::render_site(quiet = TRUE)
rmarkdown::render_site("index.Rmd")
rmarkdown::render_site("experiment.Rmd")
#rmarkdown::render_site("seqQC.Rmd")
rmarkdown::render_site("assembly.Rmd")
rmarkdown::render_site("annotation.Rmd")
rmarkdown::render_site("tools.Rmd")
