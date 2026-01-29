#################Loading########################################################
library(GSVA)
library(GSA)
library(glue)
library(Seurat)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(AUCell)
library(pheatmap)
library(patchwork)
library(flextable)


obj <- readRDS("03b-aml-PB-integrated-old.rds")
head(obj@meta.data)
colnames(obj@meta.data)
table(obj@meta.data$timing, obj@meta.data$patient)








set.seed(123)

aml.BM <- readRDS("./05a-aml-BM-annotated-complete.rds")
aml.BM$azimuthNames <- NA
aml.BM$azimuthNames[grepl("healthy", aml.BM$sample)] <- read.table("../data/byproducts/azimuth_pred_HC_BM.tsv", header = T, sep = "\t")$predicted.celltype.l2
blasts.BM <- aml.BM[,grepl("BC", aml.BM$cellType) | aml.BM$azimuthNames == "HSC"]
blasts.BM$cellType <- as.character(blasts.BM$cellType)
blasts.BM$cellType[grepl("healthy", blasts.BM$sample)] <- "HSC"
Idents(blasts.BM) <- blasts.BM$cellType <- factor(blasts.BM$cellType, levels = c("HSC", "BC1", "BC2", "BC3", "BC4"))
blasts.BM$cluster <- paste0(blasts.BM$cellType, " ", blasts.BM$timing)
blasts.BM$cluster[grepl("HSC", blasts.BM$cluster)] <- "HSC"
blasts.BM$cluster <- factor(blasts.BM$cluster, levels = c("HSC",
                                                          "BC1 pre", "BC1 post",
                                                          "BC2 pre", "BC2 post",
                                                          "BC3 pre", "BC3 post",
                                                          "BC4 pre", "BC4 post"))


aml.PB <- readRDS("../data/byproducts/05b-aml-PB-annotated-complete.rds")
aml.PB$azimuthNames <- NA
aml.PB$azimuthNames[grepl("healthy", aml.PB$sample)] <- read.table("../data/byproducts/azimuth_pred_HC_PB.tsv", header = T, sep = "\t")$predicted.celltype.l2
blasts.PB <- aml.PB[,grepl("PC", aml.PB$cellType) | aml.PB$azimuthNames == "HSPC"]
blasts.PB$cellType <- as.character(blasts.PB$cellType)
blasts.PB$cellType[grepl("healthy", blasts.PB$sample)] <- "HSPC"
Idents(blasts.PB) <- blasts.PB$cellType <- factor(blasts.PB$cellType, levels = c("HSPC", "PC"))
blasts.PB$cluster <- paste0(blasts.PB$cellType, " ", blasts.PB$timing)
blasts.PB$cluster[grepl("HSPC", blasts.PB$cluster)] <- "HSPC"
blasts.PB$cluster <- factor(blasts.PB$cluster, levels = c("HSPC",
                                                          "PC pre", "PC post"))

#Loading gene sets from Zeng et al.'s paper:
gene.sets <- GSA.read.gmt("../data/metadata/AMLCellType_Genesets.gmt")

Zeng.quiescent <- gene.sets$genesets[gene.sets$geneset.names == "LSPC-Quiescent"][[1]]%>%stringi::stri_remove_empty()
Zeng.primed <- gene.sets$genesets[gene.sets$geneset.names == "LSPC-Primed-Top250"][[1]]%>%stringi::stri_remove_empty()
Zeng.cycling <- gene.sets$genesets[gene.sets$geneset.names == "LSPC-Cycle-Top250"][[1]]%>%stringi::stri_remove_empty()
Zeng.gmp <- gene.sets$genesets[gene.sets$geneset.names == "GMP-like-Top250"][[1]]%>%stringi::stri_remove_empty()
Zeng.pro.mono <- gene.sets$genesets[gene.sets$geneset.names == "ProMono-like-Top250"][[1]]%>%stringi::stri_remove_empty()
Zeng.mono <- gene.sets$genesets[gene.sets$geneset.names == "Mono-like-Top250"][[1]]%>%stringi::stri_remove_empty()
Zeng.cdc <- gene.sets$genesets[gene.sets$geneset.names == "cDC-like-Top250"][[1]]%>%stringi::stri_remove_empty()

#And Ng et al.'s LSC17:
lsc17 <- c('GPR56', 'AKR1C3', 'CD34', 'NGFRAP1',
           'EMP1', 'SMIM24', 'SOCS2', 'CPXM1',
           'CDK6', 'KIAA0125', 'DPYSL3', 'MMRN1',
           'LAPTM4B', 'ARHGAP22', 'NYNRIN', 'ZBTB46',
           'DNMT3B')

our9 <- c("NRIP1", "ERG", "MYB", "CD34", "KMT2A", "PROM1", "NOTCH1", "BCL11A", "BCL2L1")

blasts.BM <- ScaleData(blasts.BM, features = union(Zeng.cdc, y = c(Zeng.cycling, Zeng.gmp, Zeng.mono, Zeng.primed, Zeng.pro.mono, Zeng.quiescent, lsc17)))
ave.exp <- AverageExpression(blasts.BM, slot = "scale.data", features = union(Zeng.cdc, y = c(Zeng.cycling, Zeng.gmp, Zeng.mono, Zeng.primed, Zeng.pro.mono, Zeng.quiescent, lsc17)))

#################GSVA###########################################################

gsva.score.our9 <- gsva(blasts.BM@assays$RNA@data, 
                         list(our9),
                         method="ssgsea",abs.ranking=F)

gsva.score.our9.PB <- gsva(blasts.PB@assays$RNA@data, 
                         list(our9),
                         method="ssgsea",abs.ranking=F)

gsva.score.lsc17 <- gsva(blasts.BM@assays$RNA@data, 
                         list(lsc17),
                         method="ssgsea",abs.ranking=F)

gsva.score.quiescent <- gsva(blasts.BM@assays$RNA@data, 
                             list(Zeng.quiescent),
                             method="ssgsea",abs.ranking=F)

gsva.score.quiescent.PB <- gsva(blasts.PB@assays$RNA@data, 
                             list(Zeng.quiescent),
                             method="ssgsea",abs.ranking=F)

gsva.score.primed <- gsva(blasts.BM@assays$RNA@data, 
                          list(Zeng.primed),
                          method="ssgsea",abs.ranking=F)

gsva.score.cycling <- gsva(blasts.BM@assays$RNA@data, 
                           list(Zeng.cycling),
                           method="ssgsea",abs.ranking=F)

gsva.score.cycling.PB <- gsva(blasts.PB@assays$RNA@data, 
                           list(Zeng.cycling),
                           method="ssgsea",abs.ranking=F)

gsva.score.gmp <- gsva(blasts.BM@assays$RNA@data, 
                       list(Zeng.gmp),
                       method="ssgsea",abs.ranking=F)

gsva.score.promono <- gsva(blasts.BM@assays$RNA@data, 
                           list(Zeng.pro.mono),
                           method="ssgsea",abs.ranking=F)

gsva.score.mono <- gsva(blasts.BM@assays$RNA@data, 
                        list(Zeng.mono),
                        method="ssgsea",abs.ranking=F)

gsva.score.cdc <- gsva(blasts.BM@assays$RNA@data, 
                       list(Zeng.cdc),
                       method="ssgsea",abs.ranking=F)


gsva.score <- data.frame(row.names = colnames(blasts.BM), cell = colnames(blasts.BM),
                         quiescent = gsva.score.quiescent@x,
                         primed = gsva.score.primed@x, cycling = gsva.score.cycling@x,
                         gmp = gsva.score.gmp@x, promono = gsva.score.promono@x,
                         mono = gsva.score.mono@x, cDC = gsva.score.cdc@x,
                         our9 = gsva.score.our9@x)

#Saving scores to blasts object:
blasts.BM$LSC17 <- gsva.score.lsc17@x
blasts.BM$quiescent <- gsva.score$quiescent
blasts.BM$primed <- gsva.score$primed
blasts.BM$cycling <- gsva.score$cycling
blasts.BM$gmp <- gsva.score$gmp
blasts.BM$promono <- gsva.score$promono
blasts.BM$mono <- gsva.score$mono
blasts.BM$cDC <- gsva.score$cDC
blasts.BM$our9 <- gsva.score$our9
saveRDS(blasts.BM, "../data/byproducts/05a-aml-BM-blasts-scored-complete.rds", compress = F)

#Saving scores to blasts object:
blasts.PB$quiescent <- gsva.score.quiescent.PB@x
blasts.PB$cycling <- gsva.score.cycling.PB@x
blasts.PB$our9 <- gsva.score.our9.PB@x
saveRDS(blasts.PB, "../data/byproducts/05b-aml-PB-blasts-scored-complete.rds", compress = F)

#Doing some quick t tests:
df <- blasts.BM@meta.data

#Setting up cluster comparisons for later plots:
comparisons <- list(c("BC1", "BC2"), c("BC1", "BC3"), c("BC2", "BC3"))
vs.hc <- list(c("BC1", "HSC"), c("BC2", "HSC"), c("BC3", "HSC"))
timing.comparisons <- list(c("BC1 post", "BC1 pre"),
                           c("BC2 post", "BC2 pre"),
                           c("BC3 post", "BC3 pre"))



#Plotting something that tells us how each cluster performs:
#our9
ggsave(filename = "../results/graphs/boxplot-aml-BM-our9-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "our9", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = vs.hc) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("our9 score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-our9-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df[df$response != "Healthy",], x = "cellType", y = "our9", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("our9 score"))

#LSC17
ggsave(filename = "../results/graphs/boxplot-aml-BM-LSC17-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "LSC17", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("LSC17 score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-LSC17-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "LSC17", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("LSC17 score"))

#Quiescent:
ggsave(filename = "../results/graphs/boxplot-aml-BM-quiescent-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "quiescent", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Quiescence score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-quiescent-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "quiescent", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Quiescence score"))

#Primed:
ggsave(filename = "../results/graphs/boxplot-aml-BM-primed-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "primed", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Priming score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-primed-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "primed", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Priming score"))

#Cycling:
ggsave(filename = "../results/graphs/boxplot-aml-BM-cycling-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "cycling", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Cycling score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-cycling-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "cycling", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Cycling score"))

#GMP:
ggsave(filename = "../results/graphs/boxplot-aml-BM-gmp-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "gmp", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("GMP score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-gmp-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "gmp", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("GMP score"))

#Pro-mono:
ggsave(filename = "../results/graphs/boxplot-aml-BM-promono-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "promono", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Pro-mono score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-promono-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "promono", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Pro-mono score"))

#Mono:
ggsave(filename = "../results/graphs/boxplot-aml-BM-mono-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "mono", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Mono score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-mono-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "mono", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Mono score"))

#cDC:
ggsave(filename = "../results/graphs/boxplot-aml-BM-cdc-score-cluster-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "cDC", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("cDC score"))

ggsave(filename = "../results/graphs/boxplot-aml-BM-cdc-score-cluster-response-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "cDC", outlier.shape = NA, facet.by = "response",
                        color = "cellType", bxp.errorbar = T) +
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("cDC score"))


#Comparing each cluster (pre) to HSCs:
#Quiescence:
ggsave(filename = "../results/graphs/boxplot-aml-BM-quiescent-score-cluster-pre-vs-hsc-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df[df$timing == "pre" | df$cellType == "HSC",], x = "cellType", y = "quiescent", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = vs.hc) +
         labs(color="Cluster") +
         xlab("Cluster") + ylab("Quiescence score"))
df[df$timing == "pre" | df$cellType == "HSC",] %>% wilcox_test(quiescent ~ cellType) %>%
  add_significance() %>%
  write.table("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-pre-blasts-vs-hsc-complete.tsv",
              sep = "\t", row.names = F)
#by response:
for (resp in c("responders", "non-responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-quiescent-score-cluster-", resp, "-pre-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "quiescent", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Quiescence score"))
        df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(quiescent ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-pre-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
        
        
      })
  
  try(silent = T, 
      expr = {
        
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-quiescent-score-cluster-", resp, "-post-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "quiescent", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Quiescence score"))
        df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(quiescent ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-post-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
        
      })
  
  
  
}

#Priming:
ggsave(filename = "../results/graphs/boxplot-aml-BM-priming-score-cluster-pre-vs-hsc-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df[df$timing == "pre" | df$cellType == "HSC",], x = "cellType", y = "primed", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = vs.hc) +
         labs(color="Cluster") +
         xlab("Cluster") + ylab("Priming score"))
df[df$timing == "pre" | df$cellType == "HSC",] %>% wilcox_test(primed ~ cellType) %>%
  add_significance() %>% 
  write.table(file = "../results/tables/ssgsea-wilcox-aml-BM-priming-scores-pre-blasts-vs-hsc-complete.tsv",
              sep = "\t", row.names = F)
#by response:
for (resp in c("responders", "non-responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-priming-score-cluster-", resp, "-pre-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "primed", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Priming score"))
        df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(primed ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-pre-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
        
        
        
      }
      
  )
  
  try(silent = T, 
      expr = {
        
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-priming-score-cluster-", resp, "-post-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "primed", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Priming score"))
        df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(primed ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-post-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
      })
  
  
}

#Cycling:
ggsave(filename = "../results/graphs/boxplot-aml-BM-cycling-score-cluster-pre-vs-hsc-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df[df$timing == "pre" | df$cellType == "HSC",], x = "cellType", y = "cycling", outlier.shape = NA,
                        color = "cellType", bxp.errorbar = T) + 
         stat_compare_means(comparisons = vs.hc) +
         labs(color="Cluster") +
         xlab("Cluster") + ylab("Cycling score"))
df[df$timing == "pre" | df$cellType == "HSC",] %>% wilcox_test(cycling ~ cellType) %>%
  add_significance() %>% 
  write.table(file = "../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-pre-blasts-vs-hsc-complete.tsv",
              sep = "\t", row.names = F)
#by response:
for (resp in c("responders", "non-responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-cycling-score-cluster-", resp, "-pre-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "cycling", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Cycling score"))
        df[(df$timing == "pre" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(cycling ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-pre-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
      }
  )
  
  try(silent = T,
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-cycling-score-cluster-", resp, "-post-vs-hsc-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",], 
                                x = "cellType", y = "cycling", outlier.shape = NA,
                                color = "cellType", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = vs.hc) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Cycling score"))
        df[(df$timing == "post" & df$response == resp) | df$cellType == "HSC",] %>% wilcox_test(cycling ~ cellType) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-post-blasts-", resp, "-vs-hsc-complete.tsv"),
                      sep = "\t", row.names = F)
      }
  )
  
  
  
}


df.PB <- blasts.PB@meta.data
df.PB$cluster.response <- paste0(df.PB$cluster, " ", df.PB$response)
df.PB$cluster.response <- factor(df.PB$cluster.response, levels = c("PC pre non-responders", "PC pre responders", 
                                                              "PC post non-responders", "PC post responders",
                                                              "HSPC controls"))

#Comparing post & pre for each cluster, separating by response:
#our9:
for (resp in c("non-responders", "responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-our9-score-cluster-post-vs-pre-", resp, "-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$response == resp | df$cluster == "HSC",], 
                                x = "cluster", y = "quiescent", outlier.shape = NA,
                                color = "cluster", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = timing.comparisons) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Our 9 score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        df[df$response == resp,] %>% wilcox_test(quiescent ~ cluster) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-our9-scores-post-vs-pre-", resp, "-complete.tsv"),
                      sep = "\t", row.names = F)
      })
  
  
  
}
#our9:
for (resp in c("non-responders", "responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-PB-our9-score-cluster-post-vs-pre-", resp, "-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$response == resp | df$cluster == "HSPC",], 
                                x = "cluster", y = "quiescent", outlier.shape = NA,
                                color = "cluster", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = timing.comparisons) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Our 9 score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        df[df$response == resp,] %>% wilcox_test(quiescent ~ cluster) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-PB-our9-scores-post-vs-pre-", resp, "-complete.tsv"),
                      sep = "\t", row.names = F)
      })
  
  
  
}

#Quiescence:
for (resp in c("non-responders", "responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-quiescent-score-cluster-post-vs-pre-", resp, "-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$response == resp | df$cluster == "HSC",], 
                                x = "cluster", y = "quiescent", outlier.shape = NA,
                                color = "cluster", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = timing.comparisons) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Quiescence score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        compare_means(data = df[df$response == resp,], quiescent ~ cluster)%>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-post-vs-pre-", resp, "-complete.tsv"),
                      sep = "\t", row.names = F)
      })
  
  
  
}

#Priming:
for (resp in c("non-responders", "responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-priming-score-cluster-post-vs-pre-", resp, "-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$response == resp | df$cluster == "HSC",], 
                                x = "cluster", y = "primed", outlier.shape = NA,
                                color = "cluster", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = timing.comparisons) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Priming score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        df[df$response == resp,] %>% wilcox_test(primed ~ cluster) %>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-post-vs-pre-", resp, "-complete.tsv"),
                      sep = "\t", row.names = F)
        
      })
  
  
  
}

#Cycling:
for (resp in c("non-responders", "responders")) {
  
  try(silent = T, 
      expr = {
        ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-cycling-score-cluster-post-vs-pre-", resp, "-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$response == resp | df$cluster == "HSC",], 
                                x = "cluster", y = "cycling", outlier.shape = NA,
                                color = "cluster", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = timing.comparisons) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Cycling score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        compare_means(data = df[df$response == resp,], cycling ~ cluster)%>%
          add_significance() %>%
          write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-post-vs-pre-", resp, "-complete.tsv"),
                      sep = "\t", row.names = F)
      })
  
  
  
}


#Comparing pre Refractory vs Sensitive for each cluster:
#our9:
df$cluster.response <- paste0(df$cluster, " ", df$response)
df$cluster.response <- factor(df$cluster.response, levels = c("HSC controls",
                                                              "BC1 pre non-responders", "BC1 pre responders", 
                                                              "BC2 pre non-responders", "BC2 pre responders",
                                                              "BC3 pre non-responders", "BC3 pre responders",
                                                              "BC1 post non-responders", "BC1 post responders", 
                                                              "BC2 post non-responders", "BC2 post responders",
                                                              "BC3 post non-responders", "BC3 post responders"))
ggsave(filename = paste0("../results/graphs/boxplot-aml-BM-our9-score-cluster-refractory-vs-sensitive-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df[df$timing == "pre" & !df$cluster == "HSC",], 
                                x = "cluster.response", y = "our9", outlier.shape = NA,
                                color = "cluster.response", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = list(c("BC1 pre non-responders", "BC1 pre responders"),
                                                       c("BC2 pre non-responders", "BC2 pre responders"),
                                                       c("BC3 pre non-responders", "BC3 pre responders"))) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Our 9 score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
df[df$timing == "pre" & !df$cluster == "HSC",] %>% droplevels() %>% wilcox_test(our9 ~ cluster.response) %>%
  add_significance() %>%
  write.table(paste0("../results/tables/ssgsea-wilcox-aml-BM-our9-scores-refractory-vs-sensitive-complete.tsv"),
                      sep = "\t", row.names = F)
      
ggsave(filename = paste0("../results/graphs/boxplot-aml-PB-our9-score-cluster-refractory-vs-sensitive-complete.jpeg"),
               width = 7, height = 7, unit = "in", dpi = 300,
               plot = ggboxplot(df.PB[df.PB$timing == "pre" & !df.PB$cluster == "HSPC",], 
                                x = "cluster.response", y = "our9", outlier.shape = NA,
                                color = "cluster.response", bxp.errorbar = T) + 
                 stat_compare_means(comparisons = list(c("PC pre non-responders", "PC pre responders"))) +
                 labs(color="Cluster") +
                 xlab("Cluster") + ylab("Our 9 score") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
df.PB[df.PB$timing == "pre" & !df.PB$cluster == "HSPC",] %>% droplevels() %>% wilcox_test(our9 ~ cluster.response) %>%
  add_significance() %>%
  write.table(paste0("../results/tables/ssgsea-wilcox-aml-PB-our9-scores-refractory-vs-sensitive-complete.tsv"),
                      sep = "\t", row.names = F)
      


#Saving summary statistics for clusters in tables:
##Gene set scores:
#all:
data.frame(celltype = levels(df$cellType),
                    quiescence = paste0(aggregate(x = df$quiescent, by = list(df$cellType), FUN = mean)$x %>% format(digits = 3), "(", 
                                        quiescence.sd = aggregate(x = df$quiescent, by = list(df$cellType), FUN = sd)$x %>% format(digits = 3), ")"),
                    
                    priming = paste0(aggregate(x = df$primed, by = list(df$cellType), FUN = mean)$x %>% format(digits = 3), "(",
                                     aggregate(x = df$primed, by = list(df$cellType), FUN = sd)$x %>% format(digits = 3), ")"),
                    
                    cycling = paste0(aggregate(x = df$cycling, by = list(df$cellType), FUN = mean)$x %>% format(digits = 3), "(",
                                     cycling.sd = aggregate(x = df$cycling, by = list(df$cellType), FUN = sd)$x %>% format(digits = 3), ")")) %>%
  write.table(sep = "\t", file = "../results/tables/ssgsea-aml-BM-scores-blasts-all-hsc-mean-sd-complete.tsv", row.names = F)

#separated by timing:
data.frame(cluster = levels(df$cluster),
           quiescence = paste0(aggregate(x = df$quiescent, by = list(df$cluster), FUN = mean)$x %>% format(digits = 3), "(", 
                               quiescence.sd = aggregate(x = df$quiescent, by = list(df$cluster), FUN = sd)$x %>% format(digits = 3), ")"),
           
           priming = paste0(aggregate(x = df$primed, by = list(df$cluster), FUN = mean)$x %>% format(digits = 3), "(",
                            aggregate(x = df$primed, by = list(df$cluster), FUN = sd)$x %>% format(digits = 3), ")"),
           
           cycling = paste0(aggregate(x = df$cycling, by = list(df$cluster), FUN = mean)$x %>% format(digits = 3), "(",
                            cycling.sd = aggregate(x = df$cycling, by = list(df$cluster), FUN = sd)$x %>% format(digits = 3), ")")) %>%
  write.table(sep = "\t", file = "../results/tables/ssgsea-aml-BM-scores-blasts-pre-post-hsc-mean-sd-complete.tsv", row.names = F)


df$timing.response <- paste0(df$cluster, " ", df$response) %>% as.factor()
#separated by timing, by response:
data.frame(cluster = levels(df$timing.response),
           quiescence = paste0(aggregate(x = df$quiescent, by = list(df$timing.response), FUN = mean)$x %>% format(digits = 3), "(", 
                               quiescence.sd = aggregate(x = df$quiescent, by = list(df$timing.response), FUN = sd)$x %>% format(digits = 3), ")"),
           
           priming = paste0(aggregate(x = df$primed, by = list(df$timing.response), FUN = mean)$x %>% format(digits = 3), "(",
                            aggregate(x = df$primed, by = list(df$timing.response), FUN = sd)$x %>% format(digits = 3), ")"),
           
           cycling = paste0(aggregate(x = df$cycling, by = list(df$timing.response), FUN = mean)$x %>% format(digits = 3), "(",
                            cycling.sd = aggregate(x = df$cycling, by = list(df$timing.response), FUN = sd)$x %>% format(digits = 3), ")")) %>%
  write.table(sep = "\t", file = "../results/tables/ssgsea-aml-BM-scores-blasts-timing-response-mean-sd-complete.tsv", row.names = F)


#post vs pre for each cluster:
#our9:
df[df$response == "non-responders",] %>%
  wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"),
                                 c("BC2 post", "BC2 pre"), 
                                 c("BC3 post", "BC3 pre")),
              our9 ~ cluster) %>% add_significance() %>%
  write.table(sep = "\t", row.names = F,
              file = "../results/tables/ssgsea-wilcox-aml-BM-our9-scores-refractory-post-vs-pre-complete.tsv")
df[df$response == "responders",] %>%
  wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"),
                                 c("BC2 post", "BC2 pre"), 
                                 c("BC3 post", "BC3 pre")),
              our9 ~ cluster) %>% add_significance() %>%
  write.table(sep = "\t", row.names = F,
              file = "../results/tables/ssgsea-wilcox-aml-BM-our9-scores-sensitive-post-vs-pre-complete.tsv")

#quiescence:
df %>%
  wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"),
                                 c("BC2 post", "BC2 pre"), 
                                 c("BC3 post", "BC3 pre")),
              quiescent ~ cluster) %>% add_significance() %>%
  write.table(sep = "\t", row.names = F,
              file = "../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-post-vs-pre-complete.tsv")

#priming:
df %>%
  wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"),
                                 c("BC2 post", "BC2 pre"), 
                                 c("BC3 post", "BC3 pre")),
              primed ~ cluster) %>% add_significance() %>%
  write.table(sep = "\t", row.names = F,
              file = "../results/tables/ssgsea-wilcox-aml-BM-priming-scores-post-vs-pre-complete.tsv")

#cycling:
df %>%
  wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"),
                                 c("BC2 post", "BC2 pre"), 
                                 c("BC3 post", "BC3 pre")),
              cycling ~ cluster) %>% add_significance() %>%
  write.table(sep = "\t", row.names = F,
              file = "../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-post-vs-pre-complete.tsv")



#Saving the statistical comparisons in tables:
##Quiescence score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(quiescent ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(quiescent ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-post-complete.xlsx")
#all:
df %>% wilcox_test(quiescent ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-quiescence-scores-complete.xlsx")

##Priming score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(primed ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(primed ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-post-complete.xlsx")
#all:
df %>% wilcox_test(primed ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-priming-scores-complete.xlsx")

##Cycling score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(cycling ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(cycling ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-post-complete.xlsx")
#all:
df %>% wilcox_test(cycling ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-aml-BM-cycling-scores-complete.xlsx")






#Doing kind of the same thing, but for Patient 6, comparing pre to post:
pt6 <- blasts.BM[,grepl("Patient6", blasts.BM$sample)]

#Trying out some semi-bulk stuff first:
pt6 <- ScaleData(pt6, features = union(Zeng.quiescent, y = c(Zeng.primed, Zeng.cycling)))
ave.exp <- AverageExpression(pt6, slot = "scale.data", group.by = c("ident", "timing"),
                             features = union(Zeng.quiescent, y = c(Zeng.primed, Zeng.cycling)))

pt6.scores <- gsva(ave.exp$RNA, annotation = c('Zeng.quiescent', 'Zeng.primed', 'Zeng.cycling'),
                   list(Zeng.quiescent, Zeng.primed, Zeng.cycling),
                   method="ssgsea", abs.ranking=F, parallel.sz = 20)
rownames(pt6.scores) <- c('Zeng.quiescent', 'Zeng.primed', 'Zeng.cycling')

ggsave(filename = "../results/graphs/heatmap-pt6-BM-zeng-score-cluster-timing-complete.jpeg",
       width = 5, height = 2.2, unit = "in", dpi = 300,
       plot = pheatmap(pt6.scores, cluster_cols = F, gaps_col = c(2,4,6), angle_col = 90))



#Then single cell, as done above:
gsva.score <- gsva(as.matrix(pt6@assays$RNA@data), 
                   list(Zeng.quiescent, Zeng.primed, Zeng.cycling),
                   method="ssgsea",abs.ranking=F)


rownames(gsva.score) <- c("Zeng.quiescent", "Zeng.primed", "Zeng.cycling")


pt6$quiescent <- gsva.score[rownames(gsva.score) == "Zeng.quiescent"]
pt6$primed <- gsva.score[rownames(gsva.score) == "Zeng.primed"]
pt6$cycling <- gsva.score[rownames(gsva.score) == "Zeng.cycling"]

#Saving feature plots showing score changes with time:
#Quiescence:
ggsave(filename = "../results/graphs/feature-plot-pt6-quiescent-score-timing-complete.jpeg", 
       width = 12, height = 6, units = "in", dpi = 300,
       plot = FeaturePlot(pt6, features = "quiescent", cols = c("green", "red"), 
                          split.by = "timing", shape.by = "cellType") +
         theme(legend.position = "right", legend.direction = "vertical") +
         guides(shape = guide_legend(override.aes = list(size = 4))))

#Priming:
ggsave(filename = "../results/graphs/feature-plot-pt6-primed-score-timing-complete.jpeg", 
       width = 12, height = 6, units = "in", dpi = 300,
       plot = FeaturePlot(pt6, features = "primed", cols = c("green", "red"), 
                          split.by = "timing", shape.by = "cellType")  +
         theme(legend.position = "right", legend.direction = "vertical") +
         guides(shape = guide_legend(override.aes = list(size = 4))))

#Cycling:
ggsave(filename = "../results/graphs/feature-plot-pt6-cycling-score-timing-complete.jpeg", 
       width = 12, height = 6, units = "in", dpi = 300,
       plot = FeaturePlot(pt6, features = "cycling", cols = c("green", "red"), 
                          split.by = "timing", shape.by = "cellType") +
         guides(shape = guide_legend(override.aes = list(size = 4))) +
         theme(legend.position = "right", legend.direction = "vertical"))

df <- pt6@meta.data



#Quiescent:
ggsave(filename = "../results/graphs/boxplot-pt6-BM-quiescent-score-cluster-timing-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "quiescent", outlier.shape = NA, panel.labs = list(timing = c("Pre", "Post", "HC")),
                        color = "cellType", bxp.errorbar = T, facet.by = "timing") + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Quiescence score"))

#Primed:
ggsave(filename = "../results/graphs/boxplot-pt6-BM-primed-score-cluster-timing-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "primed", outlier.shape = NA, panel.labs = list(timing = c("Pre", "Post", "HC")),
                        color = "cellType", bxp.errorbar = T, facet.by = "timing") + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Priming score"))

#Cycling:
ggsave(filename = "../results/graphs/boxplot-pt6-BM-cycling-score-cluster-timing-complete.jpeg",
       width = 7, height = 7, unit = "in", dpi = 300,
       plot = ggboxplot(df, x = "cellType", y = "cycling", outlier.shape = NA, panel.labs = list(timing = c("Pre", "Post", "HC")),
                        color = "cellType", bxp.errorbar = T, facet.by = "timing") + 
         stat_compare_means(comparisons = comparisons) +
         labs(color="Cluster") +
         xlab("Blast Cluster") + ylab("Cycling score"))



#Saving the statistical comparisons in tables:
##Quiescence score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(quiescent ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-quiescence-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(quiescent ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-quiescence-scores-post-complete.xlsx")

##Priming score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(primed ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-priming-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(primed ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-priming-scores-post-complete.xlsx")

##Cycling score:
#pre:
df[df$timing == "pre",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(cycling ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-cycling-scores-pre-complete.xlsx")
#post:
df[df$timing == "post",] %>% mutate(cellType = as.character(cellType)) %>% wilcox_test(cycling ~ cellType) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-cycling-scores-post-complete.xlsx")



#Comparing post-pre for each cluster in Patient 6:
df$clister <- paste0(df$cellType, " ", df$timing)

#Quiescence:
df[df$response == "non-responders",] %>% wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"), c("BC2 post", "BC2 pre"), c("BC3 post", "BC3 pre")), quiescent ~ clister) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  set_caption("Qiescence score comparison for each cluster - Patient 6") %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-quiescent-scores-post-vs-pre-complete.xlsx")




#Priming:
df[df$response == "non-responders",] %>% wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"), c("BC2 post", "BC2 pre"), c("BC3 post", "BC3 pre")), primed ~ clister) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  set_caption("Priming score comparison for each cluster - Patient 6") %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-priming-scores-post-vs-pre-complete.xlsx")




#Cycling:
df[df$response == "non-responders",] %>% wilcox_test(comparisons = list(c("BC1 post", "BC1 pre"), c("BC2 post", "BC2 pre"), c("BC3 post", "BC3 pre")), cycling ~ clister) %>% add_significance() %>% flextable() %>% 
  colformat_double(j = c("statistic"), digits = 0) %>%
  colformat_double(j = c("p", "p.adj"), digits = 3) %>%
  set_header_labels(values = list(p.adj.signif = "signif.")) %>%
  autofit() %>%
  delete_columns(j = 1) %>%
  set_caption("Cycling score comparison for each cluster - Patient 6") %>%
  align(j = 1:8, part = "all", align = "left") %>% 
  YesSiR::exportxlsx("../results/tables/ssgsea-wilcox-pt6-cycling-scores-post-vs-pre-complete.xlsx")
