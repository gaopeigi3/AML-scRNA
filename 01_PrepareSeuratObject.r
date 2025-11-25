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


#################Loading Data###################################################
#Loading data:
Pt1.pre <- Read10X("../data/from-cellranger/Patient1_pre/count/filtered_feature_bc_matrix/")
Pt1.post <- Read10X("../data/from-cellranger/Patient1_post/count/filtered_feature_bc_matrix/")
Pt2.pre <- Read10X("../data/from-cellranger/Patient2_pre/count/filtered_feature_bc_matrix/")
Pt2.post <- Read10X("../data/from-cellranger/Patient2_post/count/filtered_feature_bc_matrix/")
Pt3.pre <- Read10X("../data/from-cellranger/Patient3_pre/count/filtered_feature_bc_matrix/")
Pt3.post <- Read10X("../data/from-cellranger/Patient3_post/count/filtered_feature_bc_matrix/")
Pt4.pre <- Read10X("../data/from-cellranger/Patient4_pre/count/filtered_feature_bc_matrix/")
Pt4.post <- Read10X("../data/from-cellranger/Patient4_post/count/filtered_feature_bc_matrix/")
Pt5.post <- Read10X("../data/from-cellranger/Patient5_post/count/filtered_feature_bc_matrix/")
Pt6.pre.1 <- Read10X("../data/from-cellranger/Patient6_pre1/count/filtered_feature_bc_matrix/")
Pt6.pre.2 <- Read10X("../data/from-cellranger/Patient6_pre2/count/filtered_feature_bc_matrix/")
Pt6.pre.2.cd34 <- Read10X("../data/from-cellranger/Patient6_pre2_CD34/count/filtered_feature_bc_matrix/")
Pt6.post <- Read10X("../data/from-cellranger/Patient6_post/count/filtered_feature_bc_matrix/")
Pt7.pre <- Read10X("../data/from-cellranger/Patient7_pre/count/filtered_feature_bc_matrix/")
Pt8.post <- Read10X("../data/from-cellranger/Patient8_post/count/filtered_feature_bc_matrix/")
Pt9.pre <- Read10X("../data/from-cellranger/Patient9_pre/count/filtered_feature_bc_matrix/")
Pt10.pre <- Read10X("../data/from-cellranger/Patient10_pre/count/filtered_feature_bc_matrix/")
Pt10.post <- Read10X("../data/from-cellranger/Patient10_post/count/filtered_feature_bc_matrix/")
Pt11.pre <- Read10X("../data/from-cellranger/Patient11_pre/count/filtered_feature_bc_matrix/")
Pt11.post <- Read10X("../data/from-cellranger/Patient11_post/count/filtered_feature_bc_matrix/")
Pt12.pre <- Read10X("../data/from-cellranger/Patient12_pre/count/filtered_feature_bc_matrix/")
Pt13.post <- Read10X("../data/from-cellranger/Patient13_post/count/filtered_feature_bc_matrix/")
Pt14.post <- Read10X("../data/from-cellranger/Patient14_post/count/filtered_feature_bc_matrix/")
Pt15.pre <- Read10X("../data/from-cellranger/Patient15_pre/count/filtered_feature_bc_matrix/")
Pt15.post <- Read10X("../data/from-cellranger/Patient15_post/count/filtered_feature_bc_matrix/")

#Samples from June 21, 2024:
Pt16.pre <- Read10X("../data/from-cellranger/Patient16_pre/count/filtered_feature_bc_matrix/") #PB
Pt17.pre <- Read10X("../data/from-cellranger/Patient17_pre/count/sample_filtered_feature_bc_matrix/") #BM
Pt18.pre <- Read10X("../data/from-cellranger/Patient18_pre/count/sample_filtered_feature_bc_matrix/") #BM
Pt19.pre <- Read10X("../data/from-cellranger/Patient19_pre/count/sample_filtered_feature_bc_matrix/") #PB
Pt20.pre <- Read10X("../data/from-cellranger/Patient20_pre/count/sample_filtered_feature_bc_matrix/") #PB
Pt21.post <- Read10X("../data/from-cellranger/Patient21_post/count/sample_filtered_feature_bc_matrix/") #PB
Pt22.pre <- Read10X("../data/from-cellranger/Patient22_pre/count/sample_filtered_feature_bc_matrix/") #BM

#HCs:
healthy1 <- Read10X("../data/from-cellranger/healthy1/sample_feature_bc_matrix/")
healthy2 <- Read10X("../data/from-cellranger/healthy2/sample_feature_bc_matrix/")
healthy3 <- Read10X("../data/from-cellranger/healthy3/sample_feature_bc_matrix/")
BM.HC <- readRDS("../data/byproducts/nlbm-object.rds")

#Creating Seurat objects:
Pt1.pre <- CreateSeuratObject(counts = Pt1.pre, min.cells = 3, min.features = 100, project = "AML")
Pt1.post <- CreateSeuratObject(counts = Pt1.post, min.cells = 3, min.features = 100, project = "AML")
Pt2.pre <- CreateSeuratObject(counts = Pt2.pre, min.cells = 3, min.features = 100, project = "AML")
Pt2.post <- CreateSeuratObject(counts = Pt2.post, min.cells = 3, min.features = 100, project = "AML")
Pt3.pre <- CreateSeuratObject(counts = Pt3.pre, min.cells = 3, min.features = 100, project = "AML")
Pt3.post <- CreateSeuratObject(counts = Pt3.post, min.cells = 3, min.features = 100, project = "AML")
Pt4.pre <- CreateSeuratObject(counts = Pt4.pre, min.cells = 3, min.features = 100, project = "AML")
Pt4.post <- CreateSeuratObject(counts = Pt4.post, min.cells = 3, min.features = 100, project = "AML")
Pt5.post <- CreateSeuratObject(counts = Pt5.post, min.cells = 3, min.features = 100, project = "AML")
Pt6.pre.1 <- CreateSeuratObject(counts = Pt6.pre.1, min.cells = 3, min.features = 100, project = "AML")
Pt6.pre.2 <- CreateSeuratObject(counts = Pt6.pre.2, min.cells = 3, min.features = 100, project = "AML")
Pt6.pre.2.cd34 <- CreateSeuratObject(counts = Pt6.pre.2.cd34, min.cells = 3, min.features = 100, project = "AML")
Pt6.post <- CreateSeuratObject(counts = Pt6.post, min.cells = 3, min.features = 100, project = "AML")
Pt7.pre <- CreateSeuratObject(counts = Pt7.pre, min.cells = 3, min.features = 100, project = "AML")
Pt8.post <- CreateSeuratObject(counts = Pt8.post, min.cells = 3, min.features = 100, project = "AML")
Pt9.pre <- CreateSeuratObject(counts = Pt9.pre, min.cells = 3, min.features = 100, project = "AML")
Pt10.pre <- CreateSeuratObject(counts = Pt10.pre, min.cells = 3, min.features = 100, project = "AML")
Pt10.post <- CreateSeuratObject(counts = Pt10.post, min.cells = 3, min.features = 100, project = "AML")
Pt11.pre <- CreateSeuratObject(counts = Pt11.pre, min.cells = 3, min.features = 100, project = "AML")
Pt11.post <- CreateSeuratObject(counts = Pt11.post, min.cells = 3, min.features = 100, project = "AML")
Pt12.pre <- CreateSeuratObject(counts = Pt12.pre, min.cells = 3, min.features = 100, project = "AML")
Pt13.post <- CreateSeuratObject(counts = Pt13.post, min.cells = 3, min.features = 100, project = "AML")
Pt14.post <- CreateSeuratObject(counts = Pt14.post, min.cells = 3, min.features = 100, project = "AML")
Pt15.pre <- CreateSeuratObject(counts = Pt15.pre, min.cells = 3, min.features = 100, project = "AML")
Pt15.post <- CreateSeuratObject(counts = Pt15.post, min.cells = 3, min.features = 100, project = "AML")
Pt16.pre <- CreateSeuratObject(counts = Pt16.pre, min.cells = 3, min.features = 100, project = "AML")
Pt17.pre <- CreateSeuratObject(counts = Pt17.pre, min.cells = 3, min.features = 100, project = "AML")
Pt18.pre <- CreateSeuratObject(counts = Pt18.pre, min.cells = 3, min.features = 100, project = "AML")
Pt19.pre <- CreateSeuratObject(counts = Pt19.pre, min.cells = 3, min.features = 100, project = "AML")
Pt20.pre <- CreateSeuratObject(counts = Pt20.pre, min.cells = 3, min.features = 100, project = "AML")
Pt21.post <- CreateSeuratObject(counts = Pt21.post, min.cells = 3, min.features = 100, project = "AML")
Pt22.pre <- CreateSeuratObject(counts = Pt22.pre, min.cells = 3, min.features = 100, project = "AML")
healthy1 <- CreateSeuratObject(counts = healthy1, min.cells = 3, min.features = 100, project = "AML")
healthy2 <- CreateSeuratObject(counts = healthy2, min.cells = 3, min.features = 100, project = "AML")
healthy3 <- CreateSeuratObject(counts = healthy3, min.cells = 3, min.features = 100, project = "AML")
healthy4 <- BM.HC[,grepl("NLBM4", colnames(BM.HC))]
healthy5 <- BM.HC[,grepl("NLBM5", colnames(BM.HC))]
healthy6 <- BM.HC[,grepl("NLBM6", colnames(BM.HC))]

aml <- merge(x = Pt1.pre, y = c(Pt1.post, Pt2.pre, Pt2.post, Pt3.pre, Pt3.post,
                                Pt4.pre, Pt4.post, Pt5.post, Pt6.pre.1, Pt6.pre.2,
                                Pt6.pre.2.cd34, Pt6.post, Pt7.pre, Pt8.post, Pt9.pre,
                                Pt10.pre, Pt10.post, Pt11.pre, Pt11.post,
                                Pt12.pre, Pt13.post, Pt14.post, Pt15.pre, Pt15.post,
                                Pt16.pre, Pt17.pre, Pt18.pre, Pt19.pre,
                                Pt20.pre, Pt21.post, Pt22.pre,
                                healthy1, healthy2, healthy3, 
                                healthy4, healthy5, healthy6), 
               add.cell.ids = c("Patient1_pre", "Patient1_post",
                                "Patient2_pre", "Patient2_post",
                                "Patient3_pre", "Patient3_post",
                                "Patient4_pre", "Patient4_post",
                                "Patient5_post", 
                                "Patient6_pre1", "Patient6_pre2",
                                "Patient6_pre2_cd34", "Patient6_post",
                                "Patient7_pre", "Patient8_post",
                                "Patient9_pre",
                                "Patient10_pre", "Patient10_post",
                                "Patient11_pre", "Patient11_post",
                                "Patient12_pre", "Patient13_post",
                                "Patient14_post",
                                "Patient15_pre", "Patient15_post",
                                "Patient16_pre", "Patient17_pre",
                                "Patient18_pre","Patient19_pre",
                                "Patient20_pre",
                                "Patient21_post","Patient22_pre",
                                "healthy1_control", "healthy2_control", "healthy3_control",
                                "healthy4_control", "healthy5_control", "healthy6_control"),
               project = "aml")

head(aml@meta.data)

#Automated alternative for loading data:
#Loading data from folders inside "from-cellranger":
# count.mat <- list()
# for (i in list.dirs("../data/from-cellranger", recursive = F, full.names = F) ) {
#   count.mat[i] <- Read10X(paste0("../data/from-cellranger/",i, "/count/filtered_feature_bc_matrix"))
# }
# 
# #Creating list of Seurat objects:
# seurat.obj <- list()
# for (j in names(count.mat)) {
#   seurat.obj[[j]] <- CreateSeuratObject(counts = count.mat[[j]], min.cells = 3, min.features = 100)
# }
# 
# #Merging list items:
# aml <- merge(x = seurat.obj[[1]], y = seurat.obj[2:length(seurat.obj)], add.cell.ids = names(seurat.obj),
#              project = "aml")

#Saving file:
saveRDS(aml, "../data/byproducts/00-aml-raw-complete.rds", compress = F)

rm(list = setdiff(ls(), "aml"))
gc()


#################QC#############################################################
aml <-  readRDS("../data/byproducts/00-aml-raw-complete.rds")

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
aml$source[aml$sample == "Patient3_post" | aml$sample == "Patient4_pre" | aml$sample == "Patient6_pre1" | aml$sample == "Patient6_pre2_cd34" | aml$sample == "Patient7_pre" | aml$sample == "Patient9_pre" | aml$sample == "Patient15_pre" | aml$sample == "Patient15_post" | aml$sample == "Patient16_pre" | aml$sample == "Patient19_pre" | aml$sample == "Patient20_pre" | aml$sample == "Patient21_post" | aml$sample == "healthy1_control" | aml$sample == "healthy2_control" | aml$sample == "healthy3_control"] <- "PB"

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
counts <- GetAssayData(aml, slot = "counts")

#Summing up non-zeros and returning genes with 10 or more values:
keep <- rowSums(counts) >= 10

#Creating final(?) Seurat object:
aml <- CreateSeuratObject(counts = counts[keep,],
                            meta.data = aml@meta.data, project = "AML")
head(aml@meta.data)
aml$nCount_RNA <- NULL
aml$nFeature_RNA <- NULL
aml$orig.ident <- NULL
remove(list = setdiff(ls(), "aml"))
gc()

#Saving filtered file:
saveRDS(aml, "../data/byproducts/01-aml-filtered-complete.rds", compress = F)

