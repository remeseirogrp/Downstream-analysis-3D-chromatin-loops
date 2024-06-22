##CRAN packages

install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)
install.packages("tidyverse")
library(tidyverse)
install.packages("stringr")
library(stringr)
install.packages("ggpubr")
library(ggpubr)
install.packages("rstatix")
library(rstatix)
install.packages("circlize")
library(circlize)
install.packages("devtools")
library(devtools)

#Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("GenomicRanges")
library(GenomicRanges)
BiocManager::install("rtracklayer")
library(rtracklayer)
BiocManager::install("BSgenome")
library(BSgenome)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("InteractionSet")
library(InteractionSet)
BiocManager::install("GenomicInteractions")
library(GenomicInteractions)
BiocManager::install("edgeR")
library(edgeR)
BiocManager::install("limma")
library(limma)
BiocManager::install("DESeq2")
library(DESeq2)
install.packages("RColorBrewer")
library(RColorBrewer)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

#devtools packages
install_github("PhanstielLab/Sushi")
library(Sushi)
install_github("aryeelab/diffloop")
library(diffloop)














