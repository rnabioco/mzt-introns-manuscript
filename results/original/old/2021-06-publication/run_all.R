# run all markdowns
library(here)
library(rmarkdown)
base_dir <- here("results/2021-06-publication")
to_run <- dir(base_dir, pattern = ".Rmd")
walk(to_run, render)