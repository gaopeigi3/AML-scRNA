#################Environment####################################################
#Loading libraries:
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()


#################Clustering#####################################################

#For BM samples:
aml.BM <- readRDS("../data/byproducts/03a-aml-BM-integrated-complete.rds")

#Graph-based clustering:
aml.BM <- FindNeighbors(aml.BM, dims = 1:40, 
                        verbose = T, reduction = "harmony")
aml.BM <- FindClusters(aml.BM, resolution = seq(0,1,0.1), random.seed = 123, 
                       verbose = T)

#Non-linear dimensional reduction:
aml.BM <- RunUMAP(aml.BM, dims = 1:40, reduction = "harmony", seed.use = 123,
                  return.model = T, repulsion.strength = 2,
                  verbose = T)
aml.BM <- RunTSNE(aml.BM, dims = 1:40, verbose = T) 

Idents(aml.BM) <- aml.BM$RNA_snn_res.0.4

ggsave(filename = "../results/graphs/umap-aml-BM-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))
ggsave(filename = "../results/graphs/tsne-aml-BM-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "tsne", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))

gc()

saveRDS(aml.BM, "../data/byproducts/04a-aml-BM-clustered-complete.rds", compress = F)

remove(list = ls())

#For PBMC samples:
aml.PB <- readRDS("../data/byproducts/03b-aml-PB-integrated-complete.rds")

#Graph-based clustering:
aml.PB <- FindNeighbors(aml.PB, dims = 1:40, 
                        verbose = T, reduction = "harmony")
aml.PB <- FindClusters(aml.PB, resolution = seq(0,1,0.1), random.seed=123, 
                       verbose = T)

#Non-linear dimensional reduction:
aml.PB <- RunUMAP(aml.PB, dims = 1:40, reduction = "harmony", seed.use = 123,
                  return.model = T, repulsion.strength = 2,
                  verbose = T)
aml.PB <- RunTSNE(aml.PB, dims = 1:40, verbose = T) #Here, I can specify features = "CD24".


ggsave(filename = "../results/graphs/umap-aml-PB-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "umap", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))
ggsave(filename = "../results/graphs/tsne-aml-PB-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "tsne", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))

gc()

saveRDS(aml.PB, "../data/byproducts/04b-aml-PB-clustered-complete.rds", compress = F)
