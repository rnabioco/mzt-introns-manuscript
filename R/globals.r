# globals shared across markdown docs

library(readr)
library(valr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(purrr)
library(R.utils)
library(viridis)
library(jsonlite)
library(rtracklayer)
library(ComplexHeatmap)
library(tximport) 
library(DESeq2)
library(tibble)
library(here)
library(ggrepel)
library(readxl)
library(eisaR)
library(doParallel)


#### Paths ####

project_dir <- here()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

#### functions ####
source(file.path(project_dir, "R", "utils.r"))

ggplot2::theme_set(theme_cowplot())


#### Annotations ####

annots <- list(
  drosophila = list(
    gtf  = "~/Projects/shared_dbases/annotation/drosophila/Drosophila_melanogaster.BDGP6.84.gtf",
    genome = "~/Projects/shared_dbases/genomes/drosophila/Drosophila_melanogaster.BDGP6.dna.toplevel.fa",
    chroms = "~/Projects/shared_dbases/genomes/drosophila/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.fai"
  )
)
