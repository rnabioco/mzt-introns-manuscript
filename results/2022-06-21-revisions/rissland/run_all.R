# run all markdowns
library(here)
library(rmarkdown)
base_dir <- here("results/2022-06-21-revisions")
setwd(file.path(base_dir, "other-species"))
to_run <- c("01_coral.Rmd",
            "01_zebrafish.Rmd",
            "01_xenopus.Rmd")

walk(to_run, render, output_dir = "html")

# setwd(file.path(base_dir, "rissland"))
# to_run <- c("00.Rmd",
#   "01_fly.Rmd",
#   "02.Rmd",
#   "03.Rmd",
#   "nanopore.Rmd",
#   "gene-plots.Rmd")
# walk(to_run, render)