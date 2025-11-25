#################Loading########################################################
Sys.getenv("CONDA_DEFAULT_ENV")

library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)
library(SeuratData)
library(SeuratDisk)
# library(Azimuth)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)
#Enable parallelization
plan("multisession", workers = 8)
#Clearing memory:
gc()


# target_dir <- "../nas-storage2/GEO_AML/cellranger/Pt4_pre_/"
# files <- list.files(target_dir, full.names = TRUE)
# rename_files <- function(file) {
#   if (grepl("barcodes", file)) {
#     file.rename(file, file.path(target_dir, "barcodes.tsv.gz"))
#   } else if (grepl("features", file)) {
#     file.rename(file, file.path(target_dir, "features.tsv.gz"))
#   } else if (grepl("matrix", file)) {
#     file.rename(file, file.path(target_dir, "matrix.mtx.gz"))
#   }
# }

# lapply(files, rename_files)



# library(stringr)
# root_dir <- "../nas-storage2/GEO_AML/cellranger/"
# sample_dirs <- list.dirs(root_dir, recursive = FALSE)

# rename_to_standard <- function(sample_dir) {
#   files <- list.files(sample_dir, full.names = TRUE)
  
#   for (file in files) {
#     if (str_detect(file, "barcodes.*\\.tsv\\.gz$")) {
#       file.rename(file, file.path(sample_dir, "barcodes.tsv.gz"))
#     } else if (str_detect(file, "features.*\\.tsv\\.gz$")) {
#       file.rename(file, file.path(sample_dir, "features.tsv.gz"))
#     } else if (str_detect(file, "matrix.*\\.mtx\\.gz$")) {
#       file.rename(file, file.path(sample_dir, "matrix.mtx.gz"))
#     }
#   }
# }

# lapply(sample_dirs, rename_to_standard)



#################Loading Data###################################################
#Loading data:
library(stringr)
library(pbapply)



base_dir <- "../nas-storage2/GEO_AML/cellranger"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)


sample_ids <- basename(sample_dirs) |> str_remove("_$") 

seurat_list <- pblapply(seq_along(sample_dirs), function(i) {
  dir <- sample_dirs[i]
  sample_id <- sample_ids[i]
  
  data <- Read10X(file.path(dir))
  seu <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 100, project = "AML")
  return(seu)
})
names(seurat_list) <- sample_ids
aml <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "AML"
)


head(aml@meta.data)

# Saving file:
saveRDS(aml, "../nas-storage2/GEO_AML/00-aml-raw.rds", compress = F)

rm(list = setdiff(ls(), "aml"))
gc()

aml <- readRDS("../nas-storage2/GEO_AML/00-aml-raw.rds")


#Adding the number of genese per UMI:
aml$log10GenesPerUMI <- log10(aml$nFeature_RNA)/log10(aml$nCount_RNA)

#Cell-level filtering#
#Adding mitochondrial gene percentages:
aml$mitoRatio <- PercentageFeatureSet(aml, pattern = "^MT-")/100
plot1 <- FeatureScatter(aml, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(aml, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Renaming columns:
aml@meta.data <- aml@meta.data %>%
  rename(seq.folder = orig.ident,
         nUMI = nCount_RNA,
         nGene = nFeature_RNA,
         prject = orig.ident)
head(aml@meta.data)

#Creating a sample column:
aml$sample <- sub("(.*)_{1}(.*)", "\\1", colnames(aml))
aml$sample <- gsub("_NLBM.*", "", aml$sample)

#Creating a source column:
aml$source <- "BM"
aml$source[aml$sample  %in% c(
  "Pt1_pre", "Pt2_pre", "Pt3_pre", "Pt7_post", "Pt8_pre",
  "Pt10_pre", "Pt9_pre", "Pt17_pre1", "Pt18_pre", "Pt18_post", "Pt22_post"
)
] <- "PB"


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
aml@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 200) +
  geom_hline(yintercept = 2000)

head(aml@meta.data)

#Sub-setting based on gene and UMI counts and mitochondrial genes:
aml <- subset(aml, subset = (nUMI >= 500) & (nGene > 200) &  (nGene < 3000) &
                  (log10GenesPerUMI > 0.78) & (mitoRatio < 0.15))


#Removing genes with zero counts:

aml.join <- JoinLayers(aml)
counts <- GetAssayData(aml.join, assay = 'RNA', layer= "counts")  
keep <- rowSums(counts) >= 10

aml.join <- CreateSeuratObject(counts = counts[keep,],
                            meta.data = aml@meta.data, project = "AML")

head(aml.join@meta.data)
aml.join$nCount_RNA <- NULL
aml.join$nFeature_RNA <- NULL
aml.join$orig.ident <- NULL

remove(list = setdiff(ls(), "aml.join"))
gc()

saveRDS(aml.join, "../nas-storage2/GEO_AML/01-aml-filtered.rds", compress = F)



#Creating final(?) Seurat object
counts <- GetAssayData(aml, slot = "counts")
keep <- rowSums(counts) >= 10
aml <- CreateSeuratObject(counts = counts[keep,],
                            meta.data = aml@meta.data, project = "AML")
head(aml@meta.data)
# aml$nCount_RNA <- NULL
# aml$nFeature_RNA <- NULL
# aml$orig.ident <- NULL
aml@meta.data$nCount_RNA <- NULL
aml@meta.data$nFeature_RNA <- NULL
aml@meta.data$orig.ident <- NULL


#Saving filtered file:
saveRDS(aml, "../nas-storage2/GEO_AML/01-aml-filtered.rds", compress = F)




#################SummarizeData########################################################


aml <- readRDS("../nas-storage2/GEO_AML/01-aml-filtered.rds")


#################Cell types#####################################################
for (i in levels(aml$sample)) {
  cells <- table(Idents(aml[,aml$sample == i])) %>% as.data.frame() 
  colnames(cells) <- c("Cell type", "Count")
  cells <- arrange(cells, desc(Count))
  cells$percentage <- paste0("(", format(round(cells$Count / sum(cells$Count) * 100, 1), nsmall = 1),"%)")
  cells$Count <- paste0(cells$Count, " ", cells$percentage)
  cells$percentage <- NULL
  write.table(cells, file = paste0("../results/tables/cell-counts-", i, ".tsv"), 
              sep = "\t", row.names = F)
}


#################Metadata#######################################################
write.table(aml@meta.data, file = "../nas-storage2/GEO_AML/metadata/metadata.tsv",
            row.names = F, sep = "\t")


# #Per sample:
# aml@meta.data |> 
#   select(patient, timing, source, timing, response, blast, category) |>
#   unique() |> 
#   write.table("../nas-storage2/GEO_AML/metadata/patient-metadata.tsv",
#               row.names = F, sep = "\t")


# #################Cluster occupancy##############################################
# cells <- list()
# for (i in levels(aml$integrated_snn_res.2.5)) {
#   cells[[i]] <- table(aml$patient[aml$integrated_snn_res.2.5 == i]) %>% as.data.frame()
#   colnames(cells[[i]]) <- c("Patient", i)
#   cells[[i]]$percentage <- paste0("(", format(round(cells[[i]][,2] / sum(cells[[i]][,2]) * 100, 1), nsmall = 1),"%)")
#   cells[[i]][,2] <- paste0(cells[[i]][,2], " ", cells[[i]]$percentage)
#   cells[[i]]$percentage <- NULL
# }
# cells <- Reduce(function(...) merge(..., all=T), cells)
# write.table(cells, file = "../results/tables/cluster-occupancy.tsv", 
#             sep = "\t", row.names = F)






#################Split-Scaling & Azimuth########################################

aml <- readRDS("../nas-storage2/GEO_AML/01-aml-filtered.rds")

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
  saveRDS(split.aml[[i]], file = paste0("../nas-storage2/GEO_AML/", names(split.aml[i]), ".rds"), compress = F)
  names(split.aml[i])
}

#Importing Azimuth names from tsv files:
for (i in names(split.aml)) {
  split.aml[[i]]$azimuthNames <- read.table(paste0("../nas-storage2/GEO_AML/", names(split.aml[i]), "_azimuth_pred.tsv"), sep = "\t", header = T)$predicted.celltype.l2
}


saveRDS(split.aml, "../nas-storage2/GEO_AML/02-aml-split-scaled.rds", compress = F)

remove(list = setdiff(ls(), "split.aml"))
gc()







#################HarmonizeData####################################################

split.aml <- readRDS("../nas-storage2/GEO_AML/02-aml-split-scaled.rds")


#################Integration####################################################

# #For PBMC sampleS:
# aml.PB <- split.aml[c("Patient3_post", "Patient6_pre1", "Patient6_pre2_cd34",
#                               "Patient9_pre", "Patient15_pre", "Patient15_post", 
#                               "healthy1_control", "healthy2_control", "healthy3_control")]

# aml.PB <- merge(aml.PB[[1]], aml.PB[2:length(aml.PB)], merge.data = T)

# DefaultAssay(aml.PB) <- "RNA"
# #Normalizing the data:
# aml.PB <- NormalizeData(aml.PB)

# #Cell-cycle scoring:
# aml.PB <- CellCycleScoring(aml.PB, s.features = cc.genes$s.genes, 
#                            g2m.features = cc.genes$g2m.genes)

# #Scaling the data:
# aml.PB <- FindVariableFeatures(aml.PB)
# aml.PB <- ScaleData(aml.PB, vars.to.regress = c("S.Score", "G2M.Score"))

# #Linear dimensional reduction:
# aml.PB <- RunPCA(aml.PB, npcs = 40, features = VariableFeatures(object = aml.PB), 
#                  verbose = T)

# #Viewing the dimplot, grouping by cell-cycle phase:
# DimPlot(aml.PB, reduction = "pca", group.by = "Phase")

# #Integrating with Harmony:
# aml.PB <- RunHarmony(aml.PB, group.by.vars = "sample", max.iter.cluster=500,lambda = 1,
#                      verbose = T, plot_convergence = T)

# saveRDS(aml.PB, compress = F,
#         "../data/byproducts/03b-aml-PB-integrated.rds")


#For BM sampleS:
bm_samples <- aml@meta.data %>%
  filter(source == "BM") %>%
  pull(sample) %>%
  unique()
aml.BM <- split.aml[c(bm_samples)]

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
library(harmony)
aml.BM <- RunHarmony(aml.BM, group.by.vars = "sample", max.iter.cluster=500,lambda = 1,
                     verbose = T, plot_convergence = T)

saveRDS(aml.BM, compress = F,
        "../nas-storage2/GEO_AML/03a-aml-BM-integrated.rds")








#################CLusterData####################################################
#Loading libraries:
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Enable parallelization
plan("multisession", workers = 24)

#Clearing memory:
gc()


#################Clustering#####################################################

#For BM samples:
aml.BM <- readRDS("../nas-storage2/GEO_AML/03a-aml-BM-integrated.rds")

#Graph-based clustering:
aml.BM <- FindNeighbors(aml.BM, dims = 1:40, 
                        verbose = T, reduction = "harmony")
aml.BM <- FindClusters(aml.BM, resolution = seq(1,2,0.2), random.seed = 123, 
                       verbose = T)

#Non-linear dimensional reduction:
aml.BM <- RunUMAP(aml.BM, dims = 1:40, reduction = "harmony", seed.use = 123,
                  return.model = T, repulsion.strength = 2,
                  verbose = T)
aml.BM <- RunTSNE(aml.BM, dims = 1:40, verbose = T) 

Idents(aml.BM) <- aml.BM$RNA_snn_res.2

ggsave(filename = "./results/graphs/umap-aml-BMc.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))
ggsave(filename = "./results/graphs/tsne-aml-BMc.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "tsne", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1))

gc()

saveRDS(aml.BM, "../nas-storage2/GEO_AML/04a-aml-BM-clustered.rds", compress = F)

remove(list = ls())

# #For PBMC samples:
# aml.PB <- readRDS("../nas-storage2/GEO_AML/03b-aml-PB-integrated.rds")

# #Graph-based clustering:
# aml.PB <- FindNeighbors(aml.PB, dims = 1:40, 
#                         verbose = T, reduction = "harmony")
# aml.PB <- FindClusters(aml.PB, resolution = seq(0,1,0.1), random.seed=123, 
#                        verbose = T)

# #Non-linear dimensional reduction:
# aml.PB <- RunUMAP(aml.PB, dims = 1:40, reduction = "harmony", seed.use = 123,
#                   return.model = T, repulsion.strength = 2,
#                   verbose = T)
# aml.PB <- RunTSNE(aml.PB, dims = 1:40, verbose = T) #Here, I can specify features = "CD24".


# ggsave(filename = "../results/graphs/umap-aml-PB.jpeg", width = 15, height = 10, dpi = 300,
#        plot = DimPlot(aml.PB, reduction = "umap", label = T, repel = T, raster = F, 
#                       label.color = "blue", pt.size = 1))
# ggsave(filename = "../results/graphs/tsne-aml-PB.jpeg", width = 15, height = 10, dpi = 300,
#        plot = DimPlot(aml.PB, reduction = "tsne", label = T, repel = T, raster = F, 
#                       label.color = "blue", pt.size = 1))

# gc()

# saveRDS(aml.PB, "../data/byproducts/04b-aml-PB-clustered.rds", compress = F)











#################AnnotationCluster####################################################
#Loading libraries:
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Enable parallelization
plan("multisession", workers = 24)

#Clearing memory:
gc()

aml.BM <- readRDS("../nas-storage2/GEO_AML/04a-aml-BM-clustered.rds")
# aml.PB <- readRDS("../data/byproducts/04b-aml-PB-clustered.rds")


#################Cluster annotation#############################################

#Bone marrow samples:
Idents(aml.BM) <- aml.BM$RNA_snn_res.2
ClusterMarkersBM <- FindAllMarkers(aml.BM, verbose = T,
                                   logfc.threshold = 0.25, return.thresh = 0.05) |>
  group_by(cluster) |> 
  arrange(desc(avg_log2FC))

write.table(ClusterMarkersBM, sep = "\t", 
            file = "./results/tables/cluster-markers-BM-2.tsv", col.names = NA)

ClusterMarkersBM <- read.table("./results/tables/cluster-markers-BM-2.tsv", header = T)

ggsave(filename = "./results/graphs/umap-aml-BMc-clustered-2.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.color = "blue", pt.size = 1) + NoLegend())
ggsave(filename = "./results/graphs/umap-aml-BM-clustered-0.4-by-sample.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", group.by = "sample", raster = F, pt.size = 1))



#Potential AML cells/blasts: 
ggsave(filename = "./results/graphs/dotplot-aml-BMc-blasts-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CD34", "CD36", "HLA-DR",
                                                               "CD13", "CD105", "CD71", "SSC")) + #Aanei et al. (2021, Frontiers)
         theme(axis.text.x = element_text(angle = 90)))
FeaturePlot(aml.BM, features = c("CD34", "CD36", "HLA-DR",
                                 "CD13", "CD105", "CD71", "SSC"),label = T)

aml.BM
#Renaming identities:
aml.BM <- RenameIdents(aml.BM, 
                     #   "8" = "T cells", "20" = "T cells",
                     #   "5" = "T cells", "0" = "T cells",
                     #   "19" = "T cells", "31" = "T cells",
                     #   "22" = "T cells", "7" = "T cells", 
                     #   "1" = "T cells", "28" = "T cells", 
                       
                     #   "3" = "NK cells", "21" = "NK cells", 
                     #   "11" = "NK cells", 
                       
                     #   "6" = "B cells", "25" = "pre B",
                     #   "27" = "Plasma",
                       
                     #   "15" = "Late Erythroid", "10" = "Late Erythroid",
                     #   "4" = "Late Erythroid", "29" = "Prog Mk",
                       
                     #   "2" = "CD14 Mono", "14" = "CD14 Mono",
                     #   "23" = "CD16 Mono", "30" = "CD14 Mono", #Cluster 30 is actually very suspicious
                     #   "12" = "DC",
                     #   "26" = "Macrophage", "17" = "GMP", 
                     #   "18" = "BC4", "16" = "BC3", 
                       "13" = "BC2", 
                       "24" = "BC1")


bc1 <- subset(aml.BM, idents = c("BC2"))
bm_samples <- aml.BM@meta.data %>%
  filter(source == "BM") %>%
  pull(sample) %>%
  unique()
cells_target <- rownames(bc1@meta.data[bc1@meta.data$sample %in% bm_samples, ])


metadata <- aml.BM@meta.data[cells_target, ]
write.csv(as.data.frame(metadata), file = "metadata.csv")
expr_matrix <- GetAssayData(bc1, assay = "RNA", slot = "counts")[, cells_target]
write.csv(as.data.frame(expr_matrix), file = "expr_matrix.csv")


ggsave(filename = "./results/graphs/umap-aml-BM-annotated-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-BM-no-hc-annotated-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM[,!grepl("healthy", aml.BM$sample)],
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-BM-hc-annotated-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM[,grepl("healthy", aml.BM$sample)], 
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

#Making cell identity tables:
#All:
nCellsBM <- table(Idents(aml.BM[,!grepl("healthy", aml.BM$sample)])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,!grepl("healthy", aml.BM$sample)]) * 100, digits = 2)
nCellsBM$CellsHC <-  as.data.frame(table(Idents(aml.BM[,grepl("healthy", aml.BM$sample)])))$Freq
nCellsBM$PercHC <- round(nCellsBM$CellsHC / sum(nCellsBM$CellsHC) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents.tsv"))
#Responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-responders.tsv"))
#Responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-responders-pre.tsv"))
#Responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-responders-post.tsv"))
#Non-responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-non-responders.tsv"))
#Non-responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-non-responders-pre.tsv"))
#Non-responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-non-responders-post.tsv"))

#Per sample: 
for (sample in levels(factor(aml.BM$sample))) {
  object <- aml.BM[,aml.BM$sample == sample]
  nCells <- table(Idents(object)) |> as.data.frame() 
  nCells <- rename(nCells, Cells = Var1, n = Freq)
  nCells$totalPerc <- round(nCells$n / ncol(object) * 100, digits = 2)
  
  #Making barplots:
  cols <- c("T cells" = "#F8766D", "NK cells" = "#E68613",
            "B cells" = "#CD9600", "pre B" = "#ABA300",
            "Plasma" = "#7CAE00", "Late Erythroid" = "#0CB702",
            "Prog Mk" = "#00BE67", "CD14 Mono" = "#00C19A",
            "CD16 Mono" = "#00BFC4", "DC" = "#00B8E7",
            "Macrophage" = "#00A9FF", "GMP" = "#8494FF", 
            "BC4" = "#C77CFF", "BC2" = "#ED68ED",
            "BC3" = "#FF61CC", "BC1" = "#FF68A1")
  
  #Patients:
  nCells <- mutate(nCells, CellType = "Cell Type") 
  ggsave(filename = paste0("../results/graphs/barplot-aml-BM-", sample, ".jpeg"),
         width = 5, height = 10, dpi = 300,
         plot = ggplot(data = nCells, aes(x = CellType, y = totalPerc, fill = Cells)) +
           geom_bar(stat = "identity") + 
           scale_fill_manual(values = cols) +
           theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))
  
  write.table(nCells, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-", sample, ".tsv"))
  ggsave(filename = paste0("../results/graphs/umap-aml-BM-annotated-1-", sample, ".jpeg"), width = 15, height = 10, dpi = 300,
         plot = DimPlot(object, reduction = "umap", label = T, repel = T, raster = F,
                        label.color = "blue", pt.size = 1) + NoLegend())
}

#Making barplots:
cols <- c("T cells" = "#F8766D", "NK cells" = "#E68613",
          "B cells" = "#CD9600", "pre B" = "#ABA300",
          "Plasma" = "#7CAE00", "Late Erythroid" = "#0CB702",
          "Prog Mk" = "#00BE67", "CD14 Mono" = "#00C19A",
          "CD16 Mono" = "#00BFC4", "DC" = "#00B8E7",
          "Macrophage" = "#00A9FF", "GMP" = "#8494FF", 
          "BC4" = "#C77CFF", "BC2" = "#ED68ED",
          "BC3" = "#FF61CC", "BC1" = "#FF68A1")

#Patients:
nCellsBM <- mutate(nCellsBM, CellType = "Cell Type") 
ggsave(filename = "../results/graphs/barplot-aml-BM-patients.jpeg",
       width = 5, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM, aes(x = CellType, y = totalPerc, fill = Cells)) +
         geom_bar(stat = "identity") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))
#Separated on X axis:
nCellsBM <- table(Idents(aml.BM), aml.BM$sample) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, Sample = Var2, n = Freq)
nCellsBM$Sample <- factor(nCellsBM$Sample, levels = c("healthy4_control", "healthy5_control", "healthy6_control",
                                                      "Patient1_pre", "Patient1_post", #Paired responders
                                                      "Patient2_pre", "Patient2_post", 
                                                      "Patient10_pre",  "Patient10_post",
                                                      "Patient11_pre", "Patient11_post",
                                                      "Patient6_pre2", "Patient6_post", #Paired non-responder
                                                      "Patient3_pre", "Patient12_pre", #Pre, non-paired responders
                                                      #Pre, non-paired non-responders
                                                      "Patient4_post", "Patient13_post", #Post, non-paired responders 
                                                      "Patient5_post", "Patient8_post", "Patient14_post")) #Post, non-paired non-responders

#HC:
ggsave(filename = "../results/graphs/barplot-aml-BM-hc-separated.jpeg",
       width = 6, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("healthy", nCellsBM$Sample),], aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pairs:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-paired-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient1_|Patient2_|Patient10|Patient11|Patient6)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pre non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-pre-non-paired-responders-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient3|Patient12)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-post-non-paired-responders-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient4|Patient13)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired non-responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-post-non-paired-non-responders-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient5|Patient8|Patient14)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))


aml.BM$cellType <- Idents(aml.BM)

saveRDS(aml.BM, "../data/byproducts/05a-aml-BM-annotated.rds", compress = F)


#Peripheral blood samples:
#Taking annotations from theg individual Seurat objects:
hc.pb <- read_rds("../data/byproducts/hc-pb-annotated.rds")
pt3.post <- read_rds("../data/byproducts/pt3-post-pb-annotated.rds")
pt4.pre <- read_rds("../data/byproducts/pt4-pre-pb-annotated.rds")
pt6.pre1 <- read_rds("../data/byproducts/pt6-pre1-pb-annotated.rds")
pt6.pre2.cd34 <- read_rds("../data/byproducts/pt6-pre2-cd34-pb-annotated.rds")
pt7.pre <- read_rds("../data/byproducts/pt7-pre-pb-annotated.rds")
pt9.pre <- read_rds("../data/byproducts/pt9-pre-pb-annotated.rds")
pt15.pre <- read_rds("../data/byproducts/pt15-pre-pb-annotated.rds")
pt15.post <- read_rds("../data/byproducts/pt15-post-pb-annotated.rds")

blud <- merge(x = hc.pb, y = c(pt3.post[,!grepl("healthy", pt3.post$sample)],
                               pt4.pre[,!grepl("healthy", pt4.pre$sample)], 
                               pt6.pre1[,!grepl("healthy", pt6.pre1$sample)],
                               pt6.pre2.cd34[,!grepl("healthy", pt6.pre2.cd34$sample)],
                               pt7.pre[,!grepl("healthy", pt7.pre$sample)],
                               pt9.pre[,!grepl("healthy", pt9.pre$sample)], 
                               pt15.pre[,!grepl("healthy", pt15.pre$sample)],
                               pt15.post[,!grepl("healthy", pt15.post$sample)]),
              merge.data = T, merge.dr = T)
aml.PB$cellType <- blud$cellType
Idents(aml.PB) <- aml.PB$cellType 

DefaultAssay(aml.PB) <- "RNA"
Idents(aml.PB) <- aml.PB$RNA_snn_res.0.6
ClusterMarkersPB <- FindAllMarkers(aml.PB, verbose = T,
                                  logfc.threshold = 0.25, return.thresh = 0.05) |>
  group_by(cluster) |> 
  arrange(desc(avg_log2FC))

write.table(ClusterMarkersPB, sep = "\t", row.names = F,
            file = "../results/tables/cluster-markers-PB-0.6.tsv")

ClusterMarkersPB <- read.table("../results/tables/cluster-markers-PB-0.6.tsv", header = T)

ggsave(filename = "../results/graphs/umap-aml-PB-clustered-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "umap", label = T, repel = T, raster = F, 
                      label.size = 8, pt.size = 1) + NoLegend())
ggsave(filename = "../results/graphs/umap-aml-PB-clustered-1-by-sample.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "umap", group.by = "sample", raster = F, pt.size = 1))

nCellsPB <- table(Idents(aml.PB), aml.PB$sample) %>%
  as.data.frame() %>% 
  mutate(cluster = Var1, sample = Var2, n = Freq, .keep = "unused") %>%
  spread(key = sample, value = n) %>%
  arrange(as.numeric(cluster))
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-clustered-0.6-idents.tsv"))
nCellsPBHC$PercHC <- round(nCellsPBHC$CellsHC / sum(nCellsPBHC$CellsHC) * 100, digits = 2)
nCellsPB <- left_join(x = nCellsPB, y = nCellsPBHC, "Cells")


#Broad cluster annotation
#T cells: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-t-cells.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("CD3D", "CD3E", "CD3G", "CD8A",
                                                               "CD8B", "CD4")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#CD4 Proliferating: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-cd4-proliferating.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("MKI67", "TOP2A", "PCLAF", "CENPF",
                                                                 "TYMS", "NUSAP1", "ASPM", "PTTG1",
                                                                 "TPX2", "RRM2")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#NK: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-nk.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("TYROBP", "GNLY", "FCER1G", "PRF1",
                                                               "CD247", "KLRF1", "CST7", "GZMB", #Azimuth
                                                               "NKG7", "KLRC1", "KLRD1", "KLRB1", 
                                                               "KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR3DL2",
                                                               "KIR3DL3", "NCAM1", "FCGR3A")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#B cells: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-b-cells.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("CD79A", "RALGPS2", "CD79B", "MS4A1", 
                                                                 "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1",
                                                                 "IGHM", "MEF2C", #Azimuth
                                                               "CD19", "B220")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#Plasma: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-plasma.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("IGHA2", "MZB1", "TNFRSF17", "DERL3",
                                                                 "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5",
                                                                 "HRASLS2", "NT5DC2")) + #Azimuth 
                                                                 
         theme(axis.text.x = element_text(angle = 90)))

#Erythroid cells: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-erythroid.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("HBD", "HBM", "AHSP", "ALAS2", 
                                                               "CA1", "SLC4A1", "IFIT1B", "TRIM58",
                                                               "SELENBP1", "TMCC2")) + #Azimuth
                                                                 
         theme(axis.text.x = element_text(angle = 90)))

#HSPCs: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-hspc.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("SPINK2", "PRSS57", "CYTL1", "EGFL7",
                                                               "GATA2", "CD34", "SMIM24", "AVP",
                                                               "MYB", "LAPTM4B")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#CD14 Mono: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-cd14-mono.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("S100A9", "CTSS", "S100A8", "LYZ",
                                                               "VCAN", "S100A12", "IL1B", "CD14",
                                                               "G0S2", "FCN1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#CD16 Mono: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-cd16-mono.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("CDKN1C", "FCGR3A", "PTPRC", "LST1",
                                                               "IER5", "MS4A7", "RHOC", "IFITM3",
                                                               "AIF1", "HES4")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Platelets:  
ggsave(filename = "../results/graphs/dotplot-aml-PB-platelets.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("PPBP", "PF4", "NRGN", "GNG11",
                                                                 "CAVIN2", "TUBB1", "CLU", "HIST1H2AC",
                                                                 "RGS18", "GP9")) + #Azimuth
                                                               
         theme(axis.text.x = element_text(angle = 90)))

#DCs: 
ggsave(filename = "../results/graphs/dotplot-aml-PB-dc.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.PB, cluster.idents = T, features = c("CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
                                                               "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3",
                                                               "HLA-DQB1", "HLA-DRB1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))


  
#Renaming identities:
aml.PB <- RenameIdents(aml.PB, "3" = "T cells",  "5" = "T cells", 
                       "1" = "T cells", "4" = "T cells", 
                       "7" = "T cells", 
                       
                       "0" = "NK cells", "22" = "NK cells",
                       "14" = "NK cells", "8" = "NK cells", 
                       "15" = "NK cells", 
                       
                       "9" = "B cells", "11" = "B cells", 
                       "23" = "B cells",
                       
                       "13" = "Erythroid",
                       "16" = "Platelets", 
                       
                       "2" = "CD14 Mono", "20" = "CD14 Mono", 
                       "21" = "CD14 Mono", 
                       
                       "12" = "CD16 Mono", "18" = "CD16 Mono",
                       
                       "17" = "DC", "19" = "DC",
                       
                       "10" = "PC1", 
                       
                       "6" = "PC2"
                       )
cols <- c("T cells" = "#F8766D", "NK cells" = "#E58700",
          "Erythroid" = "#A3A500", "B cells" = "#6BB100",
          "DC" = "#00BA38", "CD14 Mono" = "#00BF7D",
          "CD16 Mono" = "#00C0AF", "Platelets" = "#00BCD8",
          "PC2" = "#00B0F6", "PC1" = "#619CFF", 
          "PC3" = "#B983FF", "PC4" = "#E76BF3",
          "AML m (s/o) 2" = "#FF67A4")

ggsave(filename = "../results/graphs/umap-aml-PB-annotated-0.6.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "umap", label = T, repel = T, raster = F,
        label.size = 8, pt.size = 1, cols = cols) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-PB-no-hc-annotated-0.6.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB[,!grepl("healthy", aml.PB$sample)],
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 8, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-PB-hc-annotated-0.6.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB[,grepl("healthy", aml.PB$sample)], 
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 8, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-PB-hc-azimuth-0.6.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB[,grepl("healthy", aml.PB$sample)], 
                      reduction = "umap", label = T, repel = T, raster = F, group.by = 'azimuthNames', 
                      label.size = 8, pt.size = 1) + NoLegend())


aml.PB$response <- "healthy"
aml.PB$response[aml.PB$sample %in% c("Patient6_pre1", "Patient6_pre2_cd34", "Patient15_pre", "Patient15_post" )] <- "non-responders"
aml.PB$response[aml.PB$sample %in% c("Patient4_pre", "Patient7_pre", "Patient9_pre", "Patient3_post")] <- "responders"

aml.PB$timing <- sub("(.*)_(.*)", "\\2", aml.PB$sample)
aml.PB$timing[grepl("pre", aml.PB$sample)] <- "pre"
aml.PB$timing <- factor(aml.PB$timing, levels = c("pre", "post"))

#Making cell identity tables:
write.table(table(aml.PB$cellType, aml.PB$sample), sep = "\t", col.names = NA,
            file = "../results/cell-counts-pb.tsv")
#All:
nCellsPB <- table(Idents(aml.PB[,!grepl("healthy", aml.PB$sample)])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,!grepl("healthy", aml.PB$sample)]) * 100, digits = 2)
nCellsPBHC <-  as.data.frame(table(Idents(aml.PB[,grepl("healthy", aml.PB$sample)])))
nCellsPBHC <- rename(nCellsPBHC, Cells = Var1, CellsHC = Freq)
nCellsPBHC$PercHC <- round(nCellsPBHC$CellsHC / sum(nCellsPBHC$CellsHC) * 100, digits = 2)
nCellsPB <- left_join(x = nCellsPB, y = nCellsPBHC, "Cells")
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-1-idents.tsv"))
nCellsPBHC$cellType <- "Cell Type"
ggsave(filename = "../results/graphs/barplot-aml-PB-hc-all.jpeg",
       width = 6, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPBHC, aes(x = cellType, y = CellsHC, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))

#Making cell identity tables:
#All:
nCellsPB <- table(Idents(aml.PB[,!grepl("healthy", aml.PB$sample)])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,!grepl("healthy", aml.PB$sample)]) * 100, digits = 2)
nCellsPB$CellsHC <-  as.data.frame(table(Idents(aml.PB[,grepl("healthy", aml.PB$sample)])))$Freq
nCellsPB$PercHC <- round(nCellsPB$CellsHC / sum(nCellsPB$CellsHC) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents.tsv"))
#Responders:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "responders"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "responders"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-responders.tsv"))
#Responders - pre:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "responders" & aml.PB$timing == "pre"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "responders" & aml.PB$timing == "pre"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-responders-pre.tsv"))
#Responders - post:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "responders" & aml.PB$timing == "post"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "responders" & aml.PB$timing == "post"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-responders-post.tsv"))
#Non-responders:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "non-responders"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "non-responders"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-non-responders.tsv"))
#Non-responders - pre:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "non-responders" & aml.PB$timing == "pre"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "non-responders" & aml.PB$timing == "pre"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-non-responders-pre.tsv"))
#Non-responders - post:
nCellsPB <- table(Idents(aml.PB[,aml.PB$response == "non-responders" & aml.PB$timing == "post"])) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, n = Freq)
nCellsPB$totalPerc <- round(nCellsPB$n / ncol(aml.PB[,aml.PB$response == "non-responders" & aml.PB$timing == "post"]) * 100, digits = 2)
write.table(nCellsPB, sep = "\t", row.names = F,
            file = paste0("../results/tables/aml-PB-annotated-0.6-idents-non-responders-post.tsv"))

#Saving some UMAPs:
ggsave(filename = "../results/graphs/umap-aml-PB-annotated-by-response.jpeg", width = 20, height = 10, dpi = 600,
       plot = DimPlot(aml.PB[,!grepl("healthy", aml.PB$sample)], reduction = "umap", 
                      label = T, repel = T, raster = F,
                      group.by = "cellType", split.by = "response",
                      label.size = 7, pt.size = 1) + NoLegend() + ggtitle("AML by response") + theme(plot.title = element_text(hjust = 0.5)))
ggsave(filename = "../results/graphs/umap-aml-PB-annotated-responders-by-timing.jpeg", width = 20, height = 10, dpi = 300,
       plot = DimPlot(aml.PB[,!grepl("healthy", aml.PB$sample) & aml.PB$response == 'responders'], reduction = "umap", 
                      label = T, repel = T, raster = F, 
                      group.by = "cellType", split.by = "timing", 
                      label.size = 7, pt.size = 1) + NoLegend() + ggtitle('Responders') + theme(plot.title = element_text(hjust = 0.5)))
ggsave(filename = "../results/graphs/umap-aml-PB-annotated-non-responders-by-timing.jpeg", width = 20, height = 10, dpi = 600,
       plot = DimPlot(aml.PB[,!grepl("healthy", aml.PB$sample) & aml.PB$response == 'non-responders'], reduction = "umap", 
                      label = T, repel = T, raster = F,
                      group.by = "cellType", split.by = "timing", 
                      label.size = 7, pt.size = 1) + NoLegend() + ggtitle('Non-responders') + theme(plot.title = element_text(hjust = 0.5)))


#Per sample:
for (sample in levels(factor(aml.PB$sample))) {
  object <- aml.PB[,aml.PB$sample == sample]
  nCells <- table(Idents(object)) |> as.data.frame()
  nCells <- rename(nCells, Cells = Var1, n = Freq)
  nCells$totalPerc <- round(nCells$n / sum(nCells$n)*100, digits = 2)
  
  #Making barplots:
  cols <- c("T cells" = "#F8766D", "NK cells" = "#E58700",
            "Erythroid" = "#A3A500", "B cells" = "#6BB100",
            "DC" = "#00BA38", "CD14 Mono" = "#00BF7D",
            "CD16 Mono" = "#00C0AF", "Platelets" = "#00BCD8",
            "PC2" = "#00B0F6", "PC1" = "#619CFF", 
            "PC3" = "#B983FF", "PC4" = "#E76BF3",
            "AML m (s/o) 2" = "#FF67A4")
  
  #Patients:
  nCells <- mutate(nCells, CellType = "Cell Type") 
  ggsave(filename = paste0("../results/graphs/barplot-aml-PB-", sample, ".jpeg"),
         width = 5, height = 10, dpi = 300,
         plot = ggplot(data = nCells, aes(x = CellType, y = totalPerc, fill = Cells)) +
           geom_bar(stat = "identity") + 
           scale_fill_manual(values = cols) +
           theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))
  
  write.table(nCells, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-PB-annotated-1-idents-", sample, ".tsv"))
  ggsave(filename = paste0("../results/graphs/umap-aml-PB-annotated-1-", sample, ".jpeg"), width = 15, height = 10, dpi = 300,
       plot = DimPlot(object, reduction = "umap", label = T, repel = T, raster = F,
        label.size = 8, pt.size = 1, cols = cols) + NoLegend())
}

#Making barplots:
cols <- c("T cells" = "#F8766D", "NK cells" = "#E58700",
          "Erythroid" = "#A3A500", "B cells" = "#6BB100",
          "DC" = "#00BA38", "CD14 Mono" = "#00BF7D",
          "CD16 Mono" = "#00C0AF", "Platelets" = "#00BCD8",
          "PC2" = "#00B0F6", "PC1" = "#619CFF", 
          "PC3" = "#B983FF", "PC4" = "#E76BF3",
          "AML m (s/o) 2" = "#FF67A4")

#Patients:
nCellsPB <- mutate(nCellsPB, CellType = "Cell Type") 
ggsave(filename = "../results/graphs/barplot-aml-PB-patients.jpeg",
       width = 5, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPB, aes(x = CellType, y = totalPerc, fill = Cells)) +
         geom_bar(stat = "identity") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))

#Separated on X axis:
nCellsPB <- table(Idents(aml.PB), aml.PB$sample) |> as.data.frame()
nCellsPB <- rename(nCellsPB, Cells = Var1, Sample = Var2, n = Freq)
nCellsPB$Sample <- factor(nCellsPB$Sample, levels = c("healthy1_control", "healthy2_control", "healthy3_control",
                                                      "Patient6_pre1", "Patient6_pre2_cd34", #Paired non-responders
                                                      "Patient15_pre", "Patient15_post", 
                                                      "Patient4_pre",  "Patient7_pre", "Patient9_pre", #Pre, non-paired responders
                                                      #Pre, non-paired non-responders
                                                      "Patient3_post"#Post, non-paired responders 
                                                      )) #Post, non-paired non-responders
  
#HC:
ggsave(filename = "../results/graphs/barplot-aml-PB-hc-separated.jpeg",
       width = 6, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPB[grepl("healthy", nCellsPB$Sample),], aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))

#Pairs:
ggsave(filename = "../results/graphs/barplot-aml-PB-patients-paired-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPB[grepl("(Patient6|Patient15)", nCellsPB$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pre non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-PB-patients-pre-non-paired-responders-separated.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPB[grepl("(Patient4|Patient7|Patient9)", nCellsPB$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-PB-patients-post-non-paired-responders-separated.jpeg",
       width = 5, height = 10, dpi = 300,
       plot = ggplot(data = nCellsPB[grepl("Patient3", nCellsPB$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))



aml.PB$cellType <- Idents(aml.PB)

saveRDS(aml.PB, "../data/byproducts/05b-aml-PB-annotated.rds", compress = F)

