#################Environment####################################################
#Loading libraries:
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()


#################Split-Scaling & Azimuth########################################

aml <- readRDS("../data/byproducts/01-aml-filtered-complete.rds")

#Splitting Seurat object by sample to perform scaling on all samples:
#First splitting by sample:
split.aml <- SplitObject(aml, split.by = "sample")

#Scaling:
for (i in 1:length(split.aml)) {
  split.aml[[i]] <- NormalizeData(split.aml[[i]], verbose = T)
  split.aml[[i]] <- FindVariableFeatures(split.aml[[i]], verbose = T)
  split.aml[[i]] <- ScaleData(split.aml[[i]], verbose = T)
}

#Saving separate RDS files for each sample file:
for (i in names(split.aml)) {
  saveRDS(split.aml[[i]], file = paste0("../data/byproducts/", names(split.aml[i]), ".rds"), compress = F)
  names(split.aml[i])
}

#Importing Azimuth names from tsv files:
for (i in names(split.aml)) {
  split.aml[[i]]$azimuthNames <- read.table(paste0("../data/byproducts/", names(split.aml[i]), "_azimuth_pred.tsv"), sep = "\t", header = T)$predicted.celltype.l2
}

saveRDS(split.aml, "../data/byproducts/02-aml-split-scaled-complete.rds", compress = F)

remove(list = setdiff(ls(), "split.aml"))
gc()


