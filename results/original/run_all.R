# run all markdowns
library(here)
library(rmarkdown)
base_dir <- here("results/2021-09-publication")
to_run <- c("00.Rmd",
            "01_chicken.Rmd",
            "01_fly.Rmd",
            "01_xenopus.Rmd",
            "01_zebrafish.Rmd",
            "02.Rmd",
            "03.Rmd",
            "nanopore.Rmd",
            "ribo.Rmd",
            "gene-plots.Rmd")
walk(to_run, render)