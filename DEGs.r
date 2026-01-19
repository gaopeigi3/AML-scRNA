################Environment####################################################
#Loading libraries:
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(hdf5r)
library(dplyr)
library(SeuratDisk)
library(SingleCellExperiment)
# library(tidyverse)
library(patchwork)
library(Matrix)
# library(pals)
library(cowplot)
library(zellkonverter)

setwd("./AML-scRNA")

hsc <- readH5AD("adata_HC_HSC.h5ad")
assayNames(hsc)
HSC<- as.Seurat(
  sce,
  counts = NULL,
  data = "X"
)

lsc_pre_nres <- readH5AD("adata_pre_nres_HSCLSC.h5ad")
LSC_pre_nres <- as.Seurat(
  lsc_pre_nres,
  counts = NULL,
  data = "X"
)

lsc_pre_res <- readH5AD("adata_pre_res_HSCLSC.h5ad")
LSC_pre_res <- as.Seurat(
  lsc_pre_res,
  counts = NULL,
  data = "X"
)

#cell level
gene <- "DNM1L"
expr_LSC_pre_res_cell <- FetchData(
  LSC_pre_res,
  vars = gene
) %>%
  mutate(group = "LSC_pre_res")

expr_LSC_pre_nres_cell <- FetchData(
  LSC_pre_nres,
  vars = gene
) %>%
  mutate(group = "LSC_pre_nres")

expr_HSC_cell <- FetchData(
  HSC,
  vars = gene
) %>%
  mutate(group = "HSC")

compare_two_groups <- function(df1, df2, gene, g1_name, g2_name,
                               test = c("wilcox", "t")) {
  test <- match.arg(test)

  df <- rbind(df1, df2)

  mean_g1 <- mean(df[[gene]][df$group == g1_name])
  mean_g2 <- mean(df[[gene]][df$group == g2_name])

  FC <- exp(mean_g1 - mean_g2)

  if (test == "wilcox") {
    p <- wilcox.test(
      df[[gene]] ~ df$group,
      exact = FALSE
    )$p.value
  } else {
    p <- t.test(
      df[[gene]] ~ df$group,
      var.equal = FALSE
    )$p.value
  }

  list(
    FC = FC,
    p_value = p,
    n_cells = table(df$group)
  )
}

res_pre_res_vs_HSC <- compare_two_groups(
  expr_LSC_pre_res_cell,
  expr_HSC_cell,
  gene,
  g1_name = "LSC_pre_res",
  g2_name = "HSC",
  test = "wilcox"
)
res_pre_res_vs_HSC


res_pre_nres_vs_HSC <- compare_two_groups(
  expr_LSC_pre_nres_cell,
  expr_HSC_cell,
  gene,
  g1_name = "LSC_pre_nres",
  g2_name = "HSC",
  test = "wilcox"
)
res_pre_nres_vs_HSC

res_pre_res_vs_pre_nres <- compare_two_groups(
  expr_LSC_pre_res_cell,
  expr_LSC_pre_nres_cell,
  gene,
  g1_name = "LSC_pre_res",
  g2_name = "LSC_pre_nres",
  test = "wilcox"
)
res_pre_res_vs_pre_nres






expr_cell <- rbind(expr_LSC_pre_res_cell, expr_HSC_cell)
mean_LSC_pre_res_cell <- mean(expr_cell[[gene]][expr_cell$group == "LSC"])
mean_HSC_cell <- mean(expr_cell[[gene]][expr_cell$group == "HSC"])

FC_cell <- exp(mean_LSC_pre_res_cell - mean_HSC_cell)
FC_cell
p_cell <- wilcox.test(
  expr_cell[[gene]] ~ expr_cell$group,
  exact = FALSE
)$p.value
p_cell
table(expr_cell$group)


p_cell <- t.test(
  expr_cell[[gene]] ~ expr_cell$group,
  exact = FALSE
)$p.value
p_cell 






str(HSC@meta.data)
colnames(aml.BM@meta.data)
str(aml.BM@meta.data)


LSC <- subset(
  aml.BM,
  idents = c("BC1")
)
table(LSC_pre_res@meta.data$cellType)

hc_patients <- c("HCBM1", "HCBM2", "HCBM3")
aml.HC <- subset(aml.BM, subset = (patient %in% hc_patients))
HSC <- subset(
  aml.HC,
  idents = c("BC1")
)
table(HSC@meta.data$cellType)

gene <- "DNM1L"  # 或 "Drp1"，以你对象里为准
"DNM1L" %in% rownames(aml.BM)
# donor level
expr_LSC <- FetchData(
  LSC,
  vars = c(gene, "patient")
) %>%
  mutate(group = "LSC")

expr_HSC <- FetchData(
  HSC,
  vars = c(gene, "patient")
) %>%
  mutate(group = "HSC")
library(dplyr)

pb <- bind_rows(expr_LSC, expr_HSC) %>%
  group_by(patient, group) %>%
  summarise(
    mean_expr = mean(.data[[gene]]),
    .groups = "drop"
  )

mean_LSC <- mean(pb$mean_expr[pb$group == "LSC"])
mean_HSC <- mean(pb$mean_expr[pb$group == "HSC"])

FC <- exp(mean_LSC - mean_HSC)
FC
mean_LSC
mean_HSC 
wilcox.test(
  mean_expr ~ group,
  data = pb
)
FC
wilcox.test(mean_expr ~ group, data = pb)$p.value
table(pb$group)

#cell level
gene <- "DNM1L"

expr_LSC_cell <- FetchData(
  LSC,
  vars = gene
) %>%
  mutate(group = "LSC")

expr_HSC_cell <- FetchData(
  HSC,
  vars = gene
) %>%
  mutate(group = "HSC")
expr_cell <- rbind(expr_LSC_cell, expr_HSC_cell)
mean_LSC_cell <- mean(expr_cell[[gene]][expr_cell$group == "LSC"])
mean_HSC_cell <- mean(expr_cell[[gene]][expr_cell$group == "HSC"])

FC_cell <- exp(mean_LSC_cell - mean_HSC_cell)
FC_cell
p_cell <- wilcox.test(
  expr_cell[[gene]] ~ expr_cell$group,
  exact = FALSE
)$p.value
p_cell
table(expr_cell$group)


p_cell <- t.test(
  expr_cell[[gene]] ~ expr_cell$group,
  exact = FALSE
)$p.value
p_cell 