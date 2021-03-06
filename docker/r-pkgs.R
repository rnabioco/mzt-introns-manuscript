install.packages(c("BiocManager", "remotes"))

BiocManager::install(c("Biostrings",
                       "Rsamtools",
                       "GenomicRanges",
                       "GenomicAlignments", 
                       "GenomicFeatures",
                       "Gviz",
                       "eisaR", 
                       "ComplexHeatmap",
                       "DESeq2",
                       "DEXSeq",
                       "rtracklayer", 
                       "tximport", 
                       "Rsubread",
                       "biomaRt",
                       "BiocParallel"),
                     suppressUpdates=TRUE, 
                     ask=FALSE, 
                     version = "3.15")

install.packages(c('dplyr', 
                   'purrr', 
                   'readr', 
                   'stringr', 
                   'valr', 
                   'purrr', 
                   'Matrix', 
                   'RColorBrewer', 
                   'cowplot', 
                   'doParallel', 
                   'ggplot2', 
                   'ggrepel', 
                   'here', 
                   'jsonlite', 
                   'readxl', 
                   'tibble', 
                   'tidyr', 
                   'viridis',
                   'eulerr',
                   'openxlsx',
                   'ggpubr',
                   'gprofiler2',
                   'markdown',
                   'rmarkdown',
                   'sessioninfo')) 

remotes::install_github("kriemo/kentr")


