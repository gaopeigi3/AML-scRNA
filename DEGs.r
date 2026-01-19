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
expr_LSC_all_cell <- rbind(
  expr_LSC_pre_res_cell,
  expr_LSC_pre_nres_cell
) %>%
  mutate(group = "LSC")



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

res_LSC_vs_HSC <- compare_two_groups(
  expr_LSC_all_cell,
  expr_HSC_cell,
  gene,
  g1_name = "LSC",
  g2_name = "HSC",
  test = "wilcox"
)
res_LSC_vs_HSC



