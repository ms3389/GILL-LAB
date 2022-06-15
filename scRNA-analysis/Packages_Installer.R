install.packages("hdf5r")
install.packages("ggpubr")
install.packages("cowplot")
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
install.packages("Seurat")
# BiocManager::install("infercnv")
install.packages("Matrix")
BiocManager::install("DESeq2")

install.packages("devtools")

install.packages('glue')
BiocManager::install(c("monocle"))

# Next install a few more dependencies
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))

sparseMatrixStats

BiocManager::install(c('sparseMatrixStats'))
BiocManager::install(c('matrixStats'))

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)


devtools::install_github("cole-trapnell-lab/garnett")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages('R.utils')
devtools::install_github("daskelly/earlycross")

devtools::install_version("hdf5r", version = "1.0.0", repos = "http://cran.us.r-project.org")

install.packages('kableExtra')
library(kableExtra)

# remotes::install_github("cole-trapnell-lab/garnett")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# remotes::install_github("daskelly/earlycross")

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(garnett)
library(org.Hs.eg.db)
library(DESeq2)
library(DoubletFinder)




