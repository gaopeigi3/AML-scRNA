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

aml.BM <- readRDS("../data/byproducts/04a-aml-BM-clustered-complete.rds")
aml.PB <- readRDS("../data/byproducts/04b-aml-PB-clustered-complete.rds")


cols <- c("T cells" = "#F8766D", "NK cells" = "#E68613",
          "B cells" = "#CD9600", "pre B" = "#ABA300",
          "Plasma" = "#7CAE00", "Late Erythroid" = "#0CB702",
          "Prog Mk" = "#00BE67", "CD14 Mono" = "#00C19A",
          "CD16 Mono" = "#00BFC4", "DC" = "#00B8E7",
          "Macrophage" = "#00A9FF", "GMP" = "#8494FF", 
          "BC1" = "#9F8C76", "BC2" = "brown",
          "BC3" = "#5A5A5A", "BC4" = "#000050")

cols.PB <- c("T cells" = "#F8766D", "NK cells" = "#E58700",
             "Erythroid" = "#A3A500", "B cells" = "#6BB100",
             "DC" = "#00BA38", "CD14 Mono" = "#00BF7D",
             "CD16 Mono" = "#00C0AF", "Platelets" = "#00BCD8",
             "PC" = "#00ddff")

colo <- c("BC1" = "red", "BC2" = "magenta", "BC3" = "green", "BC4" = "blue",
          "PC" = "lightgrey")



#Taking annotations from the individual Seurat objects:
#Loading annotated files:

#Bone marrow samples:
hc.bm <- read_rds("../data/byproducts/annotated-samples/hc-bm-annotated-complete.rds")
pt1.pre <- read_rds("../data/byproducts/annotated-samples/pt1-pre-bm-annotated-complete.rds")
pt1.post <- read_rds("../data/byproducts/annotated-samples/pt1-post-bm-annotated-complete.rds")
pt2.pre <- read_rds("../data/byproducts/annotated-samples/pt2-pre-bm-annotated-complete.rds")
pt2.post <- read_rds("../data/byproducts/annotated-samples/pt2-post-bm-annotated-complete.rds")
pt3.pre <- read_rds("../data/byproducts/annotated-samples/pt3-pre-bm-annotated-complete.rds")
pt4.post <- read_rds("../data/byproducts/annotated-samples/pt4-post-bm-annotated-complete.rds")
pt5.post <- read_rds("../data/byproducts/annotated-samples/pt5-post-bm-annotated-complete.rds")
pt6.pre2 <- read_rds("../data/byproducts/annotated-samples/pt6-pre2-bm-annotated-complete.rds")
pt6.post <- read_rds("../data/byproducts/annotated-samples/pt6-post-bm-annotated-complete.rds")
pt8.post <- read_rds("../data/byproducts/annotated-samples/pt8-post-bm-annotated-complete.rds")
pt10.pre <- read_rds("../data/byproducts/annotated-samples/pt10-pre-bm-annotated-complete.rds")
pt10.post <- read_rds("../data/byproducts/annotated-samples/pt10-post-bm-annotated-complete.rds")
pt11.pre <- read_rds("../data/byproducts/annotated-samples/pt11-pre-bm-annotated-complete.rds")
pt11.post <- read_rds("../data/byproducts/annotated-samples/pt11-post-bm-annotated-complete.rds")
pt12.pre <- read_rds("../data/byproducts/annotated-samples/pt12-pre-bm-annotated-complete.rds")
pt13.post <- read_rds("../data/byproducts/annotated-samples/pt13-post-bm-annotated-complete.rds")
pt14.post <- read_rds("../data/byproducts/annotated-samples/pt14-post-bm-annotated-complete.rds")
pt17.pre <- read_rds("../data/byproducts/annotated-samples/pt17-pre-bm-annotated-complete.rds")
pt18.pre <- read_rds("../data/byproducts/annotated-samples/pt18-pre-bm-annotated-complete.rds")
pt22.pre <- read_rds("../data/byproducts/annotated-samples/pt22-pre-bm-annotated-complete.rds")

bone <- merge(x = hc.bm, y = c(pt1.pre[,!grepl("healthy", pt1.pre$sample)],
                               pt1.post[,!grepl("healthy", pt1.post$sample)], 
                               pt2.pre[,!grepl("healthy", pt2.pre$sample)],
                               pt2.post[,!grepl("healthy", pt2.post$sample)],
                               pt3.pre[,!grepl("healthy", pt3.pre$sample)],
                               pt4.post[,!grepl("healthy", pt4.post$sample)], 
                               pt5.post[,!grepl("healthy", pt5.post$sample)],
                               pt6.pre2[,!grepl("healthy", pt6.pre2$sample)],
                               pt6.post[,!grepl("healthy", pt6.post$sample)],
                               pt8.post[,!grepl("healthy", pt8.post$sample)],
                               pt10.pre[,!grepl("healthy", pt10.pre$sample)],
                               pt10.post[,!grepl("healthy", pt10.post$sample)],
                               pt11.pre[,!grepl("healthy", pt11.pre$sample)],
                               pt11.post[,!grepl("healthy", pt11.post$sample)],
                               pt12.pre[,!grepl("healthy", pt12.pre$sample)],
                               pt13.post[,!grepl("healthy", pt13.post$sample)],
                               pt14.post[,!grepl("healthy", pt14.post$sample)],
                               pt17.pre[,!grepl("healthy", pt17.pre$sample)],
                               pt18.pre[,!grepl("healthy", pt18.pre$sample)],
                               pt22.pre[,!grepl("healthy", pt22.pre$sample)]),
               merge.data = T, merge.dr = T)


aml.BM$cellType <- bone$cellType
Idents(aml.BM) <- aml.BM$cellType 

ggsave(filename = "../results/graphs/umap-aml-BM-annotated-1-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 8, pt.size = 1) + NoLegend())

saveRDS(aml.BM, "../data/byproducts/05a-aml-BM-annotated-complete.rds", compress = F)



#Peripheral blood samples:
hc.pb <- read_rds("../data/byproducts/annotated-samples/hc-pb-annotated.rds")
pt3.post <- read_rds("../data/byproducts/annotated-samples/pt3-post-pb-annotated-complete.rds")
pt4.pre <- read_rds("../data/byproducts/annotated-samples/pt4-pre-pb-annotated-complete.rds")
pt6.pre1 <- read_rds("../data/byproducts/annotated-samples/pt6-pre1-pb-annotated-complete.rds")
pt6.pre2.cd34 <- read_rds("../data/byproducts/annotated-samples/pt6-pre2-cd34-pb-annotated-complete.rds")
pt7.pre <- read_rds("../data/byproducts/annotated-samples/pt7-pre-pb-annotated-complete.rds")
pt9.pre <- read_rds("../data/byproducts/annotated-samples/pt9-pre-pb-annotated-complete.rds")
pt15.pre <- read_rds("../data/byproducts/annotated-samples/pt15-pre-pb-annotated-complete.rds")
pt15.post <- read_rds("../data/byproducts/annotated-samples/pt15-post-pb-annotated-complete.rds")
pt16.pre <- read_rds("../data/byproducts/annotated-samples/pt16-pre-pb-annotated-complete.rds")
pt19.pre <- read_rds("../data/byproducts/annotated-samples/pt19-pre-pb-annotated-complete.rds")
pt20.pre <- read_rds("../data/byproducts/annotated-samples/pt20-pre-pb-annotated-complete.rds")
pt21.post <- read_rds("../data/byproducts/annotated-samples/pt21-post-pb-annotated-complete.rds")

blood <- merge(x = hc.pb, y = c(pt3.post[,!grepl("healthy", pt3.post$sample)],
                               pt4.pre[,!grepl("healthy", pt4.pre$sample)], 
                               pt6.pre1[,!grepl("healthy", pt6.pre1$sample)],
                               pt6.pre2.cd34[,!grepl("healthy", pt6.pre2.cd34$sample)],
                               pt7.pre[,!grepl("healthy", pt7.pre$sample)],
                               pt9.pre[,!grepl("healthy", pt9.pre$sample)], 
                               pt15.pre[,!grepl("healthy", pt15.pre$sample)],
                               pt15.post[,!grepl("healthy", pt15.post$sample)],
                               pt16.pre[,!grepl("healthy", pt16.pre$sample)],
                               pt19.pre[,!grepl("healthy", pt19.pre$sample)],
                               pt20.pre[,!grepl("healthy", pt20.pre$sample)],
                               pt21.post[,!grepl("healthy", pt21.post$sample)]),
              merge.data = T, merge.dr = T)


aml.PB$cellType <- blood$cellType
Idents(aml.PB) <- aml.PB$cellType 


ggsave(filename = "../results/graphs/umap-aml-PB-annotated-0.3-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.PB, reduction = "umap", label = T, repel = T, raster = F,
                      label.size = 8, pt.size = 1, cols = cols.PB) + NoLegend())


saveRDS(aml.PB, "../data/byproducts/05b-aml-PB-annotated-complete.rds", compress = F)

#################Cluster annotation#############################################

#Bone marrow samples:
Idents(aml.BM) <- aml.BM$RNA_snn_res.2
ClusterMarkersBM <- FindAllMarkers(aml.BM, verbose = T,
                                   logfc.threshold = 0.25, return.thresh = 0.05) |>
  group_by(cluster) |> 
  arrange(desc(avg_log2FC))

write.table(ClusterMarkersBM, sep = "\t", 
            file = "../results/tables/cluster-markers-BM-2-complete.tsv", col.names = NA)

ClusterMarkersBM <- read.table("../results/tables/cluster-markers-BM-2-complete.tsv", header = T)

ggsave(filename = "../results/graphs/umap-aml-BM-clustered-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, 
                      label.size = 8, pt.size = 1) + NoLegend())
ggsave(filename = "../results/graphs/umap-aml-BM-clustered-2-by-sample-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", group.by = "sample", raster = F, pt.size = 1))


#Broad cluster annotation
#T cells: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-t-cells-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CD3D", "CD3E", "CD3G", "CD4",
                                                               "CD8A", "CD8B")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#NK: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-nk-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("NKG7", "KLRC1", "KLRD1", "KLRB1", 
                                                               "KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR3DL2",
                                                               "KIR3DL3", "GNLY", "NCAM1", "FCGR3A")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#B cells: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-b-cells-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("FCRL1", "FCRL2", "CD22", "ARHGAP24", #Azimuth
                                                               "BANK1", "MS4A1", "RALGPS2", "CD37",
                                                               "SWAP70", "CD79A",
                                                               "CD19", "B220")) + #Abbas Lab
         theme(axis.text.x = element_text(angle = 90)))

#pre B cells: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-pre-b-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("NPY", "LCN6", "RAG2", "HMHB1",
                                                               "ARPP21", "AKAP12", "RAG1", "C10orf10",
                                                               "CYGB", "SLC8A1-AS1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Early Erythroid cells: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-early-erythroid-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CNRIP1", "GATA2", "ITGA2B", "TFR2",
                                                               "GATA1", "KLF1", "CYTL1", "MAP7",
                                                               "FSCN1", "APOC1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Late Erythroid cells: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-late-erythroid-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CTSE", "TSPO2", "IFIT1B", "TMEM56",
                                                               "RHCE", "RHAG", "SPTA1", "ADD2",
                                                               "EPCAM", "HBG1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#CLPs: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-clp-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("ACY3", "PRSS2", "C1QTNF4", "SPINK2",
                                                               "SMIM24", "NREP", "CD34", "DNTT", "FLT3",
                                                               "SPNS3")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#GMPs: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-gmp-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("SERPINB10","RNASE3","MS4A3","PRTN3",
                                                               "ELANE","AZU1","CTSG","RNASE2","RETN",
                                                               "NPW")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Plasma: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-plasma-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("SDC1", "IGLC6", "IGLV6-57", "TNFRSF17",
                                                               "IGLV3-1", "TNFRSF13B", "IGLC7", "JSRP1",
                                                               "FCRL5", "IGKV1-5")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#DCs: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-dc-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("CLEC4C", "PROC", "SCT", "SCN9A",
                                                               "SHD", "PPM1J", "ENHO", "CLEC10A",
                                                               "LILRA4", "DNASE1L3")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Mono: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-mono-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("LYPD2", "FOLR3", "CLEC4E", "LILRA1",
                                                               "CDA", "RBP7", "CD300LF", "FPR1", 
                                                               "CD93", "MTMR11")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#CD14 Mono: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-cd14-mono-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("FOLR3", "CLEC4E", "MCEMP1", "RBP7",
                                                               "CDA", "FPR1", "CD300E", "C5AR1",
                                                               "CD93", "APOBEC3A")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#CD16 Mono: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-cd16-mono-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("LYPD2", "VMO1", "TPPP3", "C1QA",
                                                               "C5AR1", "CD300E", "GPBAR1", "LILRA1",
                                                               "HES4", "APOBEC3A")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Platelets: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-platelets-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("RGS18", "C2orf88", "TMEM40", "GP9",
                                                               "PF4", "PPBP", "DAB2", "SPARC",
                                                               "RUFY1", "F13A1")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Macrophage: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-macrophage-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DotPlot(aml.BM, cluster.idents = T, features = c("SPIC", "FABP3", "CD5L", "CCL18",
                                                               "C1QC", "C1QB", "FABP4", "C1QA",
                                                               "APOE", "SELENOP")) + #Azimuth
         theme(axis.text.x = element_text(angle = 90)))

#Potential AML cells/blasts: 
ggsave(filename = "../results/graphs/dotplot-aml-BM-blasts-complete.jpeg", width = 15, height = 10, dpi = 300,
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

ggsave(filename = "../results/graphs/umap-aml-BM-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM, reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-BM-no-hc-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
       plot = DimPlot(aml.BM[,!grepl("healthy", aml.BM$sample)],
                      reduction = "umap", label = T, repel = T, raster = F, cols = cols,
                      label.size = 7, pt.size = 1) + NoLegend())

ggsave(filename = "../results/graphs/umap-aml-BM-hc-annotated-2-complete.jpeg", width = 15, height = 10, dpi = 300,
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
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents.tsv"))
#Responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-responders.tsv"))
#Responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-responders-pre.tsv"))
#Responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-responders-post.tsv"))
#Non-responders:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-non-responders.tsv"))
#Non-responders - pre:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "pre"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-non-responders-pre.tsv"))
#Non-responders - post:
nCellsBM <- table(Idents(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"])) |> as.data.frame()
nCellsBM <- rename(nCellsBM, Cells = Var1, n = Freq)
nCellsBM$totalPerc <- round(nCellsBM$n / ncol(aml.BM[,aml.BM$response == "non-responders" & aml.BM$timing == "post"]) * 100, digits = 2)
write.table(nCellsBM, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-0.4-complete-idents-non-responders-post.tsv"))

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
  ggsave(filename = paste0("../results/graphs/barplot-aml-BM-", sample, "-0.4-complete.jpeg"),
         width = 5, height = 10, dpi = 300,
         plot = ggplot(data = nCells, aes(x = CellType, y = totalPerc, fill = Cells)) +
           geom_bar(stat = "identity") + 
           scale_fill_manual(values = cols) +
           theme_minimal(base_size = 16) + ylab("Percentage") + xlab(NULL))
  
  write.table(nCells, sep = "\t", row.names = F,
              file = paste0("../results/tables/aml-BM-annotated-1-idents-", sample, "-0.4-complete.tsv"))
  ggsave(filename = paste0("../results/graphs/umap-aml-BM-annotated-1-", sample, "-0.4-complete.jpeg"), width = 15, height = 10, dpi = 300,
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
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-0.4-complete.jpeg",
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
ggsave(filename = "../results/graphs/barplot-aml-BM-hc-separated-complete.jpeg",
       width = 6, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("healthy", nCellsBM$Sample),], aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pairs:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-paired-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient1_|Patient2_|Patient10|Patient11|Patient6)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Pre non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-pre-non-paired-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient3|Patient12)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-post-non-paired-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient4|Patient13)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))
#Post non-paired non-responders:
ggsave(filename = "../results/graphs/barplot-aml-BM-patients-post-non-paired-non-responders-separated-complete.jpeg",
       width = 10, height = 10, dpi = 300,
       plot = ggplot(data = nCellsBM[grepl("(Patient5|Patient8|Patient14)", nCellsBM$Sample),],
                     aes(x = Sample, y = n, fill = Cells)) +
         geom_bar(stat = "identity", position = "fill") + 
         scale_fill_manual(values = cols) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 90)) +
         ylab("Percentage") + xlab(NULL))



aml.BM$cellType <- Idents(aml.BM)
saveRDS(aml.BM, "../data/byproducts/05a-aml-BM-annotated-complete.rds", compress = F)



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

