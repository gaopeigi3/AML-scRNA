#################Environment####################################################
#Loading libraries:
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)
library(harmony)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 16000 * 1024^2)

#Clearing memory:
gc()

split.aml <- readRDS("../data/byproducts/02-aml-split-scaled-complete.rds")


#################Integration####################################################

#For PBMC sampleS:
aml.PB <- split.aml[c("healthy1_control", "healthy2_control", "healthy3_control",
                      "Patient15_post", "Patient15_pre", "Patient16_pre",
                      "Patient19_pre", "Patient20_pre", "Patient21_post",
                      "Patient3_post", "Patient4_pre", "Patient6_pre1",
                      "Patient6_pre2_cd34", "Patient7_pre", "Patient9_pre")]

aml.PB <- merge(aml.PB[[1]], aml.PB[2:length(aml.PB)], merge.data = T)

DefaultAssay(aml.PB) <- "RNA"
#Normalizing the data:
aml.PB <- NormalizeData(aml.PB)

#Cell-cycle scoring:
aml.PB <- CellCycleScoring(aml.PB, s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes)

#Scaling the data:
aml.PB <- FindVariableFeatures(aml.PB)
aml.PB <- ScaleData(aml.PB, vars.to.regress = c("S.Score", "G2M.Score"))

#Linear dimensional reduction:
aml.PB <- RunPCA(aml.PB, npcs = 40, features = VariableFeatures(object = aml.PB), 
                 verbose = T)

#Viewing the dimplot, grouping by cell-cycle phase:
DimPlot(aml.PB, reduction = "pca", group.by = "Phase")

#Integrating with Harmony:
aml.PB <- RunHarmony(aml.PB, group.by.vars = "sample", max.iter.cluster=500,lambda = 1,
                     verbose = T, plot_convergence = T)

saveRDS(aml.PB, compress = F,
        "../data/byproducts/03b-aml-PB-integrated-complete.rds")

#I made the following changes to the sample names:
# Removed Patient18_pre & Patient21_pre1, due to poor quality
# Changed Patient19_pre to Patient18_pre
# Changed Patient20_pre to Patient19_pre
# Changed Patient21_pre to Patient20_pre
# Changed Patient22_post to Patient21_post
# Changed Patient23_pre to Patient22_pre


#For BM sampleS:
aml.BM <- split.aml[c("healthy4_control", "healthy5_control", "healthy6_control", 
                      "Patient1_post", "Patient1_pre", "Patient10_post", "Patient10_pre",
                      "Patient11_post", "Patient11_pre", "Patient12_pre", "Patient13_post", 
                      "Patient14_post", "Patient17_pre", "Patient18_pre", "Patient2_post", 
                      "Patient2_pre", "Patient22_pre", "Patient3_pre", "Patient4_post", 
                      "Patient5_post", "Patient6_post", "Patient6_pre2", "Patient8_post")]

aml.BM <- merge(aml.BM[[1]], aml.BM[2:length(aml.BM)], merge.data = T)

DefaultAssay(aml.BM) <- "RNA"
#Normalizing the data:
aml.BM <- NormalizeData(aml.BM)

#Cell-cycle scoring:
aml.BM <- CellCycleScoring(aml.BM, s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes)

#Scaling the data:
aml.BM <- FindVariableFeatures(aml.BM)
aml.BM <- ScaleData(aml.BM, vars.to.regress = c("S.Score", "G2M.Score"))

#Linear dimensional reduction:
aml.BM <- RunPCA(aml.BM, npcs = 40, features = VariableFeatures(object = aml.BM), 
                 verbose = T)

#Viewing the dimplot, grouping by cell-cycle phase:
DimPlot(aml.BM, reduction = "pca", group.by = "Phase")

#Integrating with Harmony:
aml.BM <- RunHarmony(aml.BM, group.by.vars = "sample", max.iter.cluster=500,lambda = 1,
                     verbose = T, plot_convergence = T)

saveRDS(aml.BM, compress = F,
        "../data/byproducts/03a-aml-BM-integrated-complete.rds")
