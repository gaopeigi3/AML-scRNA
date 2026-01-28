#################Environment####################################################
#Loading libraries:
library(Seurat)
library(ggplot2)
# library(tidyverse)
library(patchwork)
library(Matrix)
# library(pals)
library(cowplot)
library(zellkonverter)
library(pheatmap)
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()
setwd('./AML-scRNA')
# aml.BM <- readRDS("./05a-aml-BM-annotated-complete.rds")

aml.BM  <- readH5AD("adata_05a_with_embeddings.h5ad")
assayNames(aml.BM )
aml.BM <- as.Seurat(
  aml.BM,
  counts = NULL,
  data = "X"
)

emb <- as.matrix(read.csv("umap_embeddings.csv", row.names = 1))
aml.BM[["umap"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(aml.BM))

rm_patients <- c("HCBM1", "HCBM2", "HCBM3")
aml.BM <- subset(aml.BM, subset = !(patient %in% rm_patients))
table(aml.BM@meta.data$patient)

aml.BM_responder <- subset(aml.BM, subset = venetoclax == "responder")
aml.BM_nonresponder <- subset(aml.BM, subset = venetoclax == "nonresponder")
aml.BM_responder_pre <- subset(aml.BM, subset = venetoclax == "responder" & treat == "pre")
aml.BM_nonresponder_pre <- subset(aml.BM, subset = venetoclax == "nonresponder" & treat == "pre")
aml.BM_responder_post <- subset(aml.BM, subset = venetoclax == "responder" & treat == "post")
aml.BM_nonresponder_post <- subset(aml.BM, subset = venetoclax == "nonresponder" & treat == "post")


colnames(aml.BM@meta.data)
str(aml.BM@meta.data)

unique(aml.BM$venetoclax)



cols <- c(
  # ---- CD4 / CD8 T cells ----
  "CD4"            = "#8dd3c7",
  "CD8 CTL"        = "#fb8072",
  "CD8 Ex"         = "#b30000",
  "CD8 Naive"      = "#fdbb84",
  "CD8 Mem"        = "#e34a33",
  "CD8 EM"         = "#f16913",

  # ---- B cells ----
  "B cells"        = "#80b1d3",
  "pre B"          = "#08306b",
  "Plasma"         = "#1c9099",

  # ---- NK ----
  "NK"             = "#fccde5",

  # ---- Progenitors / Myeloid ----
  "Prog Mk"        = "#fdb462",
  "GMP"            = "#88419d",
  "CD14 Mono"      = "#9e9ac8",
  "CD16 Mono"      = "#756bb1",
  "Macrophage"     = "#ffed6f",

  # ---- DC ----
  "DC"             = "#b3de69",

  # ---- Erythroid ----
  "Late Erythroid" = "#969696",

  # ---- AML-related BC clusters ----
  "BC1"            = "#b30000",
  "BC2"            = "#e31a1c",
  "BC3"            = "#fb9a99",
  "BC4"            = "#990000"
)





ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_pre.pdf",
  width = 15, height = 10, dpi = 600,
  plot = DimPlot(
    aml.BM_responder_pre,
    reduction = "umap",
    group.by = "dcellType",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_pre.pdf",
  width = 15, height = 10, dpi = 600,
  plot = DimPlot(
    aml.BM_nonresponder_pre,
    reduction = "umap",
    group.by = "dcellType",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_post.pdf",
  width = 15, height = 10, dpi = 600,
  plot = DimPlot(
    aml.BM_responder_post,
    reduction = "umap",
    group.by = "dcellType",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_post.pdf",
  width = 15, height = 10, dpi = 600,
  plot = DimPlot(
    aml.BM_nonresponder_post,
    reduction = "umap",
    group.by = "dcellType",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)





cell_prop_responder <- aml.BM_responder@meta.data %>%
  group_by(treat, dcellType) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) * 100)  # 百分比

ggplot(cell_prop_responder, aes(x = treat, y = freq, fill = dcellType)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols) +  # 使用你定义的颜色
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  ylab("Cell proportion (%)") +
  ggtitle("Cell type composition by treatment (Responders)") 

ggsave(
  "../results/graphs/barplot_Responders.pdf",
  width = 6, height = 4, dpi = 600
)


cell_prop_nonresponder <- aml.BM_nonresponder@meta.data %>%
  group_by(treat, dcellType) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) * 100)  # 百分比

ggplot(cell_prop_nonresponder, aes(x = treat, y = freq, fill = dcellType)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols) +  # 使用你定义的颜色
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  ylab("Cell proportion (%)") +
  ggtitle("Cell type composition by treatment (Non-responders)") 

ggsave(
  "../results/graphs/barplot_NonResponders.pdf",
  width = 6, height = 4, dpi = 600
)












# conda activate scRNA_R

library(Seurat)
library(ggplot2)
# library(tidyverse)
library(patchwork)
library(Matrix)
# library(pals)
library(cowplot)
library(zellkonverter)
library(pheatmap)
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()
setwd('./AML-scRNA')
# aml.BM <- readRDS("./05a-aml-BM-annotated-complete.rds")

aml.BM  <- readH5AD("adata_05a_with_embeddings.h5ad")

aml.BM <- as.Seurat(
  aml.BM,
  counts = NULL,
  data = "X"
)
# assayNames(aml.BM)

emb <- as.matrix(read.csv("umap_embeddings.csv", row.names = 1))
aml.BM[["umap"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(aml.BM))

rm_patients <- c("HCBM1", "HCBM2", "HCBM3")
aml.BM <- subset(aml.BM, subset = !(patient %in% rm_patients))
table(aml.BM@meta.data$patient)
table(aml.BM@meta.data$dcellType)

# 先确保是字符，不要直接在 factor 上 ifelse
aml.BM@meta.data$dcellType2 <- as.character(aml.BM@meta.data$dcellType)

aml.BM@meta.data$dcellType2[aml.BM@meta.data$dcellType2 %in% c("BC1","BC2","BC3","BC4")] <- "BC (blast clusters)"
aml.BM@meta.data$dcellType2[aml.BM@meta.data$dcellType2 %in% c("CD14 Mono","CD16 Mono")] <- "Mono"

# 最后再转回 factor（可选）
aml.BM@meta.data$dcellType2 <- factor(aml.BM@meta.data$dcellType2)

table(aml.BM@meta.data$dcellType2)



aml.BM_responder <- subset(aml.BM, subset = venetoclax == "responder")
aml.BM_nonresponder <- subset(aml.BM, subset = venetoclax == "nonresponder")
aml.BM_responder_pre <- subset(aml.BM, subset = venetoclax == "responder" & treat == "pre")
aml.BM_nonresponder_pre <- subset(aml.BM, subset = venetoclax == "nonresponder" & treat == "pre")
aml.BM_responder_post <- subset(aml.BM, subset = venetoclax == "responder" & treat == "post")
aml.BM_nonresponder_post <- subset(aml.BM, subset = venetoclax == "nonresponder" & treat == "post")


cols <- c(
  # ---- CD4 / CD8 T cells ----
  "CD4"            = "#8dd3c7",
  "CD8 CTL"        = "#fb8072",
  "CD8 Ex"         = "#b30000",
  "CD8 Naive"      = "#fdbb84",
  "CD8 Mem"        = "#e34a33",
  "CD8 EM"         = "#f16913",

  # ---- B cells ----
  "B cells"        = "#80b1d3",
  "pre B"          = "#08306b",
  "Plasma"         = "#1c9099",

  # ---- NK ----
  "NK"             = "#fccde5",

  # ---- Progenitors / Myeloid ----
  "Prog Mk"        = "#fdb462",
  "GMP"            = "#88419d",
  "Mono"      = "#9e9ac8",
  "Macrophage"     = "#ffed6f",

  # ---- DC ----
  "DC"             = "#b3de69",

  # ---- Erythroid ----
  "Late Erythroid" = "#969696",

  # ---- AML-related BC clusters ----
  "BC (blast clusters)" = "#b30000" 
)





ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_pre.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_responder_pre,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_pre_nolabel.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_responder_pre,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = FALSE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)

ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_pre.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_nonresponder_pre,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_pre_nolabel.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_nonresponder_pre,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = FALSE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)

ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_post.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_responder_post,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_responder_post_nolabel.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_responder_post,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = FALSE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)

ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_post.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_nonresponder_post,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)
ggsave(
  filename = "../results/graphs/umap-aml-BM_nonresponder_post_nolabel.pdf",
  width = 10, height = 10, dpi = 300,
  plot = DimPlot(
    aml.BM_nonresponder_post,
    reduction = "umap",
    group.by = "dcellType2",   # 👈 使用 dcellType 列
    label = FALSE,
    raster = FALSE,
    cols = cols,
    label.size = 8,
    pt.size = 1
  ) + NoLegend()
)






# ===============================
# Requirement (①): Barplot + "values for each item"
#  - produce a table with n and freq (%)
#  - save CSV
#  - optionally add text labels onto bars for editing/inspection
# ===============================
make_prop_table <- function(seu_obj, group_label) {
  df <- seu_obj@meta.data %>%
    group_by(treat, dcellType2) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(treat) %>%
    mutate(freq = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(group = group_label)

  return(df)
}

cell_prop_responder2    <- make_prop_table(aml.BM_responder,    "Responder")
cell_prop_nonresponder2 <- make_prop_table(aml.BM_nonresponder, "Non-responder")

# ---- Save values (the "each item value" they asked for)
write.csv(cell_prop_responder2,
          "../results/tables/celltype_prop_Responder_byTreat_merged.csv",
          row.names = FALSE)

write.csv(cell_prop_nonresponder2,
          "../results/tables/celltype_prop_NonResponder_byTreat_merged.csv",
          row.names = FALSE)

# ---- Barplots (with labels)
p_bar_res <- ggplot(cell_prop_responder2, aes(x = treat, y = freq, fill = dcellType2)) +
  geom_bar(stat = "identity", width = 0.8) +
  # geom_text(
  #   aes(label = sprintf("%.1f", freq)),
  #   position = position_stack(vjust = 0.5),
  #   size = 3
  # ) +
  scale_fill_manual(values = cols, drop = FALS②E) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  ylab("Cell proportion (%)") +
  ggtitle("Cell type composition by treatment (Responders; merged BC/Mono)")

ggsave(
  "../results/graphs/barplot_Responders_merged.pdf",
  width = 7, height = 4.5, dpi = 600,
  plot = p_bar_res
)

p_bar_nres <- ggplot(cell_prop_nonresponder2, aes(x = treat, y = freq, fill = dcellType2)) +
  geom_bar(stat = "identity", width = 0.8) +
  # geom_text(
  #   aes(label = sprintf("%.1f", freq)),
  #   position = position_stack(vjust = 0.5),
  #   size = 3
  # ) +
  scale_fill_manual(values = cols, drop = FALSE) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  ylab("Cell proportion (%)") +
  ggtitle("Cell type composition by treatment (Non-responders; merged BC/Mono)")

ggsave(
  "../results/graphs/barplot_NonResponders_merged.pdf",
  width = 7, height = 4.5, dpi = 600,
  plot = p_bar_nres
)

# ---- Also save a combined table (optional,便利)
cell_prop_all2 <- bind_rows(cell_prop_responder2, cell_prop_nonresponder2)
write.csv(cell_prop_all2,
          "../results/tables/celltype_prop_ALL_byTreat_merged.csv",
          row.names = FALSE)






























aml.BM_cd8 <- subset(aml.BM, subset = grepl("^CD8", dcellType))
aml.BM_cd8$group <- paste(aml.BM_cd8$venetoclax, aml.BM_cd8$treat, sep = "_")
table(aml.BM_cd8$group)
Idents(aml.BM_cd8) <- aml.BM_cd8$group



run_deg_two_group <- function(seurat_obj, group_col, ident1, ident2, logfc_cut, prefix,
                              top_n = 50, order_by = "avg_log2FC") {
  message("Running: ", ident1, " vs ", ident2, " (|log2FC|>", logfc_cut, ")")

  # ---- 仅保留两个组 ----
  obj_sub <- seurat_obj[, seurat_obj@meta.data[[group_col]] %in% c(ident1, ident2)]
  Idents(obj_sub) <- obj_sub@meta.data[[group_col]]

  # ---- 差异分析 ----
  res <- FindMarkers(
    obj_sub,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = 0,
    min.pct = 0.1
  )

  # ---- FDR & log2FC 筛选 ----
  res_filtered <- res %>%
    filter(abs(avg_log2FC) > logfc_cut & p_val_adj < 0.05)

  if (nrow(res_filtered) == 0) {
    message("⚠️ No DEGs for ", prefix, " (|log2FC|>", logfc_cut, ")")
    return(NULL)
  }

  # ---- 限制展示基因数（例如 top 50）----
  if (!is.null(top_n)) {
    res_filtered <- res_filtered %>%
      arrange(desc(abs(avg_log2FC))) %>%
      head(top_n)
  }

  genes_use <- rownames(res_filtered)

  # ---- 提取并标准化表达矩阵 ----
  mat <- GetAssayData(obj_sub, slot = "data")[genes_use, colnames(obj_sub)]
  mat_scaled <- t(scale(t(mat)))

  # ---- 注释信息 ----
  ann_col <- data.frame(Group = obj_sub@meta.data[[group_col]])
  ann_col$Group <- factor(ann_col$Group, levels = c(ident2, ident1))  # ✅ 修正处
  rownames(ann_col) <- colnames(obj_sub)

  # ---- 输出文件 ----
  out_file <- paste0("../results/heatmaps/heatmap_", prefix, "_logfc", logfc_cut, "_top", top_n, ".jpeg")

  # ---- 绘图 ----
  pheatmap(
    mat_scaled,
    annotation_col = ann_col,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,   # ✅ 显示基因名
    show_colnames = FALSE,
    border_color = NA,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    filename = out_file,
    width = 6,
    height = 10
  )

  message("✅ Saved: ", out_file, " (", nrow(res_filtered), " DEGs plotted)")
  return(res_filtered)
}




# 1) post vs pre in nonresponders
deg_1_1   <- run_deg_two_group (aml.BM_cd8, "group", "nonresponder_post", "nonresponder_pre", 1,   "1_nonres_post_vs_pre")
deg_1_1p5 <- run_deg_two_group(aml.BM_cd8, "group", "nonresponder_post", "nonresponder_pre", 1.5, "1_nonres_post_vs_pre")

# 2) post vs pre in responders
deg_2_1   <- run_deg_two_group(aml.BM_cd8, "group", "responder_post", "responder_pre", 1,   "2_res_post_vs_pre")
deg_2_1p5 <- run_deg_two_group(aml.BM_cd8, "group", "responder_post", "responder_pre", 1.5, "2_res_post_vs_pre")

# 3) post nonresponders vs post responders
deg_3_1   <- run_deg_two_group(aml.BM_cd8, "group", "nonresponder_post", "responder_post", 1,   "3_post_nonres_vs_res")
deg_3_1p5 <- run_deg_two_group(aml.BM_cd8, "group", "nonresponder_post", "responder_post", 1.5, "3_post_nonres_vs_res")

# 4) pre nonresponders vs pre responders
deg_4_1   <- run_deg_two_group(aml.BM_cd8, "group", "nonresponder_pre", "responder_pre", 1,   "4_pre_nonres_vs_res")
deg_4_1p5 <- run_deg_two_group(aml.BM_cd8, "group", "nonresponder_pre", "responder_pre", 1.5, "4_pre_nonres_vs_res")







aml.BM_cd8 <- subset(aml.BM, subset = grepl("^CD8", dcellType))
# 设置分组标识（例如 "venetoclax" + "treat" 合并成一个变量）
aml.BM_cd8$group <- paste(aml.BM_cd8$venetoclax, aml.BM_cd8$treat, sep = "_")

# 查看有哪些组合
table(aml.BM_cd8$group)

# 1) post vs pre in nonresponders
deg_1 <- FindMarkers(
  aml.BM_cd8,
  ident.1 = "nonresponder_post",
  ident.2 = "nonresponder_pre",
  group.by = "group",
  logfc.threshold = 0, # 我们自己过滤log2FC
  min.pct = 0.1
)

# 2) post vs pre in responders
deg_2 <- FindMarkers(
  aml.BM_cd8,
  ident.1 = "responder_post",
  ident.2 = "responder_pre",
  group.by = "group"
)

# 3) post nonresponders vs post responders
deg_3 <- FindMarkers(
  aml.BM_cd8,
  ident.1 = "nonresponder_post",
  ident.2 = "responder_post",
  group.by = "group"
)

# 4) pre nonresponders vs pre responders
deg_4 <- FindMarkers(
  aml.BM_cd8,
  ident.1 = "nonresponder_pre",
  ident.2 = "responder_pre",
  group.by = "group"
)

#(a) 合并结果并筛选基因
deg_list <- list(
  post_vs_pre_nonres = deg_1,
  post_vs_pre_res = deg_2,
  post_nonres_vs_res = deg_3,
  pre_nonres_vs_res = deg_4
)

# 按 log2FC 与 FDR 阈值筛选
filter_deg <- function(df, logfc_cut) {
  df %>%
    filter(abs(avg_log2FC) > logfc_cut & p_val_adj < 0.05)
}

deg_1_filtered1 <- filter_deg(deg_1, 1)
deg_1_filtered1p5 <- filter_deg(deg_1, 1.5)
deg_2_filtered1 <- filter_deg(deg_2, 1)
deg_2_filtered1p5 <- filter_deg(deg_2, 1.5)
deg_3_filtered1 <- filter_deg(deg_3, 1)
deg_3_filtered1p5 <- filter_deg(deg_3, 1.5)
deg_4_filtered1 <- filter_deg(deg_4, 1)
deg_4_filtered1p5 <- filter_deg(deg_4, 1.5)

# 提取所有显著基因并合并
all_genes <- unique(c(
  rownames(deg_1_filtered1),
  rownames(deg_2_filtered1),
  rownames(deg_3_filtered1),
  rownames(deg_4_filtered1)
))

#(c) 取表达矩阵并标准化
mat <- GetAssayData(aml.BM_cd8, slot = "data")[all_genes, ]
mat_scaled <- t(scale(t(mat)))  # 基因行标准化


library(pheatmap)

# 每个细胞的分组信息
ann_col <- data.frame(group = aml.BM_cd8$group)
rownames(ann_col) <- colnames(aml.BM_cd8)

# 颜色
ann_colors <- list(group = c(
  "responder_pre" = "#4daf4a",
  "responder_post" = "#984ea3",
  "nonresponder_pre" = "#377eb8",
  "nonresponder_post" = "#e41a1c"
))

pheatmap(
  mat_scaled,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  filename = "../results/graphs/CD8_DEG_heatmap_log2FC1.pdf", # 👈 保存路径
  width = 10,
  height = 8
)








