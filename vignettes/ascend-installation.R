## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup_bioconductor, eval = FALSE------------------------------------
#  if (!requireNamespace("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install()

## ----bioconductor_packages, eval = FALSE---------------------------------
#  bioconductor_packages <- c("Biobase", "BiocGenerics", "BiocParallel",
#                             "SingleCellExperiment", "GenomeInfoDb",
#                             "GenomeInfoDbData")
#  BiocManager::install(bioconductor_packages)

## ----install_ascend, eval = FALSE----------------------------------------
#  # Load devtools package
#  library(devtools)
#  
#  # Use devtools to install the package
#  install_github("powellgenomicslab/ascend", build_vignettes = TRUE)

## ----load_ascend---------------------------------------------------------
# Load the package in R
library(ascend)

## ----SetupNix, eval = FALSE----------------------------------------------
#  library(BiocParallel)
#  ncores <- parallel::detectCores() - 1
#  register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

## ----SetupWin, eval = FALSE----------------------------------------------
#  library(BiocParallel)
#  workers <- 3 # Number of cores on your machine - 1
#  register(SnowParam(workers = workers,
#                     type = "SOCK",
#                     progressbar = TRUE), default = TRUE)

