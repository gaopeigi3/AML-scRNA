Sys.getenv("CONDA_DEFAULT_ENV")
#################AnnotationCluster####################################################
#Loading libraries:
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)
library(SeuratDisk)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)
#Enable parallelization
plan("multisession", workers = 16)

#Clearing memory:
gc()

# aml.BM <- readRDS("../data/byproducts/04a-aml-BM-clustered-complete.rds")
aml.BM <- readRDS("../nas-storage1/AML-scRNA-seq/data/byproducts/05a-aml-BM-annotated-complete.rds")
# aml.PB <- readRDS("../data/byproducts/04b-aml-PB-clustered.rds")


head(aml.BM @meta.data)
length(unique(aml.BM@meta.data$sample))
unique(aml.BM@meta.data$sample)

# 直接用下面的seuratdisk保存会只留下前2000个基因
DefaultAssay(aml.BM) <- "RNA"  # 确保是 RNA assay
VariableFeatures(aml.BM) <- rownames(aml.BM)  # 让所有基因保留下来

SaveH5Seurat(aml.BM, filename = "aml_BM_allgenes.h5Seurat", overwrite = TRUE)
Convert("aml_BM_allgenes.h5Seurat", dest = "h5ad", overwrite = TRUE)


# 1. 导出稀疏矩阵（counts 或 normalized data）
counts <- aml.BM@assays$RNA@counts     # 原始counts矩阵（稀疏矩阵 dgCMatrix）
data <- aml.BM@assays$RNA@data         # 标准化后的表达矩阵（log normalized）

# 保存为 Matrix Market 格式（适合Python读取）
writeMM(counts, file = "aml_BM_counts.mtx")  

# 2. 保存基因名和细胞名
write.table(rownames(counts), file = "aml_BM_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts), file = "aml_BM_barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 3. 保存细胞元数据
meta.data <- aml.BM@meta.data
write.csv(meta.data, file = "aml_BM_metadata.csv", row.names = TRUE)






SaveH5Seurat(aml.BM, filename = "aml_BM.h5Seurat", overwrite = TRUE)
Convert("aml_BM.h5Seurat", dest = "h5ad", overwrite = TRUE)

length(rownames(aml.BM))  # 原始基因数
length(VariableFeatures(aml.BM))  # 将导出的基因数（默认是2000）


slotNames(aml.BM)
#  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
#  [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
#  [11] "version"      "commands"    "tools"       

aml.BM@assays

meta <- aml.BM@meta.data
str(meta)
# View(meta)
unique(meta$cellType)
levels(meta$response)  
meta$response


cols <- c("T cells" = "#F8766D", "NK cells" = "#E68613",
          "B cells" = "#CD9600", "pre B" = "#ABA300",
          "Plasma" = "#7CAE00", "Late Erythroid" = "#0CB702",
          "Prog Mk" = "#00BE67", "CD14 Mono" = "#00C19A",
          "CD16 Mono" = "#00BFC4", "DC" = "#00B8E7",
          "Macrophage" = "#00A9FF", "GMP" = "#8494FF", 
          "BC1" = "#9F8C76", "BC2" = "brown",
          "BC3" = "#5A5A5A", "BC4" = "#000050")

colo <- c("BC1" = "red", "BC2" = "magenta", "BC3" = "green", "BC4" = "blue",
          "PC" = "lightgrey")



#################Cluster annotation#############################################
Idents(aml.BM) 

#Bone marrow samples:
Idents(aml.BM) <- aml.BM$RNA_snn_res.2 #分辨度为2的结果保存为聚类簇
ClusterMarkersBM <- FindAllMarkers(aml.BM, verbose = T,
                                   logfc.threshold = 0.25, return.thresh = 0.05) |>
  group_by(cluster) |> 
  arrange(desc(avg_log2FC))

aml.BM$cellType


write.table(ClusterMarkersBM, sep = "\t", 
            file = "./results/tables/cluster-markers-BM-2-complete.tsv", col.names = NA)

ClusterMarkersBM <- read.table("./results/tables/cluster-markers-BM-2-complete.tsv", header = T)

ggsave(filename = "./results/graphs/umap-aml-BM-clustered-2-complete1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.size = 8, pt.size = 1) + NoLegend())
ggsave(filename = "./results/graphs/umap-aml-BM-clustered-2-by-sample-complete1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", group.by = "sample", raster = F, pt.size = 1))


#Potential AML cells/blasts: 
ggsave(filename = "./results/graphs/dotplot-aml-BM-blasts-complete-1.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CD34", "CD36", "HLA-DR",
                                                               "CD13", "CD105", "CD71", "SSC")) + #Aanei et al. (2021, Frontiers)
         theme(axis.text.x = element_text(angle = 90)))
FeaturePlot(aml.BM, features = c("CD34", "CD36", "HLA-DR",
                                 "CD13", "CD105", "CD71", "SSC"),label = T)

#Renaming identities:
aml.BM <- RenameIdents(aml.BM, "24" = "T cells", "1" = "T cells",
                       "11" = "T cells", "22" = "T cells",
                       "2" = "T cells", "0" = "T cells",
                       "3" = "T cells", "5" = "T cells", 
                       "7" = "T cells", "34" = "T cells", 
                       
                       "8" = "NK cells", "27" = "NK cells", 
                       "20" = "NK cells", 
                       
                       "36" = "B cells", "6" = "B cells", "28" = "B cells",
                       "26" = "pre B", "33" = "Plasma",
                       
                       "14" = "Late Erythroid", "19" = "Late Erythroid",
                       "23" = "Late Erythroid", "21" = "Late Erythroid",
                       "10" = "Late Erythroid", "12" = "Late Erythroid", 
                       "32" = "Prog Mk",
                       
                       "4" = "CD14 Mono", "29" = "CD14 Mono",
                       "35" = "CD14 Mono", "30" = "CD16 Mono", 
                       "15" = "DC",
                       "31" = "Macrophage", "18" = "GMP", 
                       "16" = "BC1", "9" = "BC2", "25" = "BC2", 
                       "13" = "BC3", "17" = "BC4")



bc1 <- subset(aml.BM, idents = c("BC1"))
DimPlot(bc1, label = TRUE)

head(bc1@meta.data)


target_samples <- c("Patient1_pre", "Patient2_pre", "Patient3_pre", 
                    "Patient6_pre", "Patient10_pre", "Patient11_pre", 
                    "Patient12_pre", "Patient17_pre", "Patient18_pre", 
                    "Patient22_pre")
cells_target <- rownames(bc1@meta.data[bc1@meta.data$sample %in% target_samples, ])

metadata <- bc1@meta.data[cells_target, ]
write.csv(metadata, file = "metadata.csv", row.names = TRUE)

expr_matrix <- GetAssayData(bc1, assay = "RNA", slot = "counts")[, cells_target]
write.csv(as.data.frame(expr_matrix), file = "expr_matrix3.csv")





# 你要的样本列表
target_samples <- c("Patient1_pre", "Patient2_pre", "Patient3_pre", 
                    "Patient6_pre", "Patient10_pre", "Patient11_pre", 
                    "Patient12_pre", "Patient17_pre", "Patient18_pre", 
                    "Patient22_pre")

# 提取所有细胞名对应的样本前缀
cell_ids <- colnames(aml.BM)  # 例如 "Patient10_pre_AAAGTC..."
sample_names <- sapply(strsplit(cell_ids, "_"), function(x) paste(x[1:2], collapse = "_"))
names(sample_names) <- cell_ids  # 把细胞名作为索引

# 获取满足两个条件的细胞：样本在目标列表，cellType 为 "BC1"
selected_cells <- names(sample_names)[
  sample_names %in% target_samples & aml.BM$cellType == "BC1"
]

# 提取 metadata 和表达数据
metadata_selected <- aml.BM@meta.data[selected_cells, ]
expression_selected <- aml.BM@assays$RNA@data[, selected_cells]

# 映射字典
replace_map <- c(
  "Patient1_pre" = "Pt11_pre",
  "Patient2_pre" = "Pt6_pre",
  "Patient3_pre" = "Pt7_pre",
  "Patient6_pre" = "Pt17_pre",
  "Patient10_pre" = "Pt12_pre",
  "Patient11_pre" = "Pt13_pre",
  "Patient12_pre" = "Pt5_pre",
  "Patient17_pre" = "Pt4_pre",
  "Patient18_pre" = "Pt15_pre",
  "Patient22_pre" = "Pt9_pre"
)
# 替换函数
replace_sample_name <- function(cell_names, map) {
  for (old in names(map)) {
    new <- map[[old]]
    cell_names <- gsub(paste0("^", old), new, cell_names)
  }
  return(cell_names)
}

# 替换 metadata 和 expression 的细胞名
rownames(metadata_selected) <- replace_sample_name(rownames(metadata_selected), replace_map)
colnames(expression_selected) <- replace_sample_name(colnames(expression_selected), replace_map)




write.csv(metadata_selected, file = "metadata5.csv", row.names = TRUE)
write.csv(as.data.frame(expression_selected), file = "expr_matrix5.csv")












#################Cluster annotation#############################################

#Bone marrow samples:
Idents(aml.BM) <- aml.BM$RNA_snn_res.2
ClusterMarkersBM <- FindAllMarkers(aml.BM, verbose = T,
                                   logfc.threshold = 0.25, return.thresh = 0.05) |>
  group_by(cluster) |> 
  arrange(desc(avg_log2FC))

write.table(ClusterMarkersBM, sep = "\t", 
            file = "./results/tables/cluster-markers-BM-2-complete.tsv", col.names = NA)

ClusterMarkersBM <- read.table("./results/tables/cluster-markers-BM-2-complete.tsv", header = T)

ggsave(filename = "./results/graphs/umap-aml-BM-clustered-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.size = 8, pt.size = 1) + NoLegend())
ggsave(filename = "./results/graphs/umap-aml-BM-clustered-2-by-sample-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", group.by = "sample", raster = F, pt.size = 1))


#Renaming identities:
aml.BM <- RenameIdents(aml.BM, "24" = "T cells", "1" = "T cells",
                       "11" = "T cells", "22" = "T cells",
                       "2" = "T cells", "0" = "T cells",
                       "3" = "T cells", "5" = "T cells", 
                       "7" = "T cells", "34" = "T cells", 
                       
                       "8" = "NK cells", "27" = "NK cells", 
                       "20" = "NK cells", 
                       
                       "36" = "B cells", "6" = "B cells", "28" = "B cells",
                       "26" = "pre B", "33" = "Plasma",
                       
                       "14" = "Late Erythroid", "19" = "Late Erythroid",
                       "23" = "Late Erythroid", "21" = "Late Erythroid",
                       "10" = "Late Erythroid", "12" = "Late Erythroid", 
                       "32" = "Prog Mk",
                       
                       "4" = "CD14 Mono", "29" = "CD14 Mono",
                       "35" = "CD14 Mono", "30" = "CD16 Mono", 
                       "15" = "DC",
                       "31" = "Macrophage", "18" = "GMP", 
                       "16" = "BC1", "9" = "BC2", "25" = "BC2", 
                       "13" = "BC3", "17" = "BC4")

ggsave(filename = "./results/graphs/umap-aml-BM-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "./results/graphs/umap-aml-BM-no-hc-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM[,!grepl("healthy", aml.BM$sample)],
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "./results/graphs/umap-aml-BM-hc-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
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
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents.tsv"))
#Responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-responders.tsv"))
#Responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-responders-pre.tsv"))
#Responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-responders-post.tsv"))
#Non-responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-non-responders.tsv"))
#Non-responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-non-responders-pre.tsv"))
#Non-responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-0.4-complete-idents-non-responders-post.tsv"))

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
  ggsave(filename = paste0("./results/graphs/barplot-aml-BM-", sample, "-0.4-complete.jpeg"),
         width = 5, height = 10, dpi = 300,
         plot = ggplot(data = nCells, aes(x = CellType, y = totalPerc, fill = Cells)) +
           geom_bar(stat = "identity") + 
           scale_fill_manual(values = cols) +
           theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))
  
  write.table(nCells, sep = "\t", row.names = F,
              file = paste0("./results/tables/aml-BM-annotated-1-idents-", sample, "-0.4-complete.tsv"))
  ggsave(filename = paste0("./results/graphs/umap-aml-BM-annotated-1-", sample, "-0.4-complete.jpeg"), width = 15, height = 10, dpi = 300,
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
ggsave(filename = "./results/graphs/barplot-aml-BM-patients-0.4-complete.jpeg",
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
ggsave(filename = "./results/graphs/barplot-aml-BM-hc-separated-complete.jpeg",
       width = 6, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("healthy", nCellsBM$Sample),], aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pairs:
ggsave(filename = "./results/graphs/barplot-aml-BM-patients-paired-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient1_|Patient2_|Patient10|Patient11|Patient6)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pre non-paired responders:
ggsave(filename = "./results/graphs/barplot-aml-BM-patients-pre-non-paired-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient3|Patient12)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired responders:
ggsave(filename = "./results/graphs/barplot-aml-BM-patients-post-non-paired-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient4|Patient13)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired non-responders:
ggsave(filename = "./results/graphs/barplot-aml-BM-patients-post-non-paired-non-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient5|Patient8|Patient14)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))



aml.BM$cellType <- Idents(aml.BM)
saveRDS(aml.BM, "./data/byproducts/05a-aml-BM-annotated-complete.rds", compress = F)

