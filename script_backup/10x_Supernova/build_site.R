#! /usr/bin/env R
#Set our working directory. 
#This helps avoid confusion if our working directory is 
#not our site because of other projects we were 
#working on at the time. 

library(rmarkdown)
setwd("report_template")

#render your sweet site. 
rmarkdown::render_site(quiet = TRUE)
