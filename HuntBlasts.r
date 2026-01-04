#################Loading########################################################
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(paletteer)
library(enrichR)

aml <- readRDS("../data/byproducts/05-aml-annotated.rds")

cols <- paletteer_d("palettesForR::Inkscape")[15:58]
names(cols) <- levels(factor(aml$azimuthNames))


#################Clustering for blasts##########################################

#Subsetting for HS(P)Cs:
hsc <- aml[,grep("HS(P)?C", aml$azimuthNames)]

#Scaling:
hsc <- ScaleData(hsc)

#Dimensionality reduction:
hsc <- RunPCA(hsc, features = VariableFeatures(object = hsc), verbose = T)

ElbowPlot(hsc, ndims = 50) #Maybe go with 40 here too...

#Graph-based clustering:
hsc <- FindNeighbors(hsc, dims = 1:40, verbose = T)
hsc <- FindClusters(hsc, verbose = T, resolution = 0.05)

#Non-linear dimensional reduction:
hsc <- RunUMAP(hsc, dims = 1:40, return.model = T, repulsion.strength = 2, verbose = T)

DimPlot(hsc, reduction = "umap", label = T, repel = T, raster = F, split.by = "timing",
        label.color = "blue", pt.size = 1)

saveRDS(hsc, "../data/byproducts/hs(p)c-clustered.rds")


#################DGE analysis between clusters##################################

hsc <- readRDS("../data/byproducts/hs(p)c-clustered.rds")

hsc.markers <- FindAllMarkers(hsc) %>% 
  filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))

hsc.markers[hsc.markers$gene %in% c("AKR1C3", "ARHGAP22", "C19orf177", "CD34", 
                                    "CDK6", "CPXM1", "DNMT3B", "DPYSL3", 
                                    "EMP1", "GPR56", "KIAA0125", "LAPTM4B",
                                    "MMRN1", "NGFRAP1", "NYNRIN", "SOCS2", "ZBTB46", #Ng's LSC17
                                    "APOC2", "CD36", "CD38", "BCL2","MEIS1", #Exploring
                                    "PPARG", "ASS1", "MYO1C", "CDC42EP4", #From our bulk RNA-seq
                                    "IL2RA", "NPDC1", "PHGDH", #Nguyen et al.'s 4-GES, along with SOCS2
                                    "TRAF3IP2", "MEF2C", "FRMD4B", "BCL11A", "COL5A1", "ELK3", "ERG", "PTK2", "SPTBN1", "ARPP-19"),] #From Vitali et al.

write.table(hsc.markers, col.names = NA, sep = "\t", file = "../results/tables/deg-clustered-hs(p)c.tsv")
write.table(hsc.markers[hsc.markers$gene %in% c("AKR1C3", "ARHGAP22", "C19orf177", "CD34", 
                                                "CDK6", "CPXM1", "DNMT3B", "DPYSL3", 
                                                "EMP1", "GPR56", "KIAA0125", "LAPTM4B",
                                                "MMRN1", "NGFRAP1", "NYNRIN", "SOCS2", "ZBTB46", #Ng's LSC17
                                                "APOC2", "CD36", "CD38", "BCL2","MEIS1", #Exploring
                                                "PPARG", "ASS1", "MYO1C", "CDC42EP4", #From our bulk RNA-seq
                                                "IL2RA", "NPDC1", "PHGDH", #Nguyen et al.'s 4-GES, along with SOCS2
                                                "TRAF3IP2", "MEF2C", "FRMD4B", "BCL11A", "COL5A1", "ELK3", "ERG", "PTK2", "SPTBN1", "ARPP-19"),],
            col.names = NA, sep = "\t", file = "../results/tables/deg-genes-clustered-hs(p)c.tsv")


table(hsc$sample, hsc$azimuthNames, hsc$integrated_snn_res.0.05)

for (i in levels(hsc$sample)) {
  ggsave(filename = paste0("../results/graphs/umap-hs(p)c-reclustered-", i, ".jpeg"),
         width = 15, height = 10, dpi = 300,
         plot = DimPlot(hsc[,hsc$sample == i], reduction = "umap", pt.size = 1) +
           xlim(-16, 7.5) + ylim(-20, 10.5) +
           labs(title = paste0(sub(x = i, pattern = "(.*)_(.*)",  "\\1"), " ",  sub(x = i, pattern = "(.*)_(.*)",  "\\2"), " (", hsc[,hsc$sample == i]$source, ")")))
}

#Trying to compare pre and post in only the suspect clusters 1/2:
for (i in levels(aml$patient)) {
  try(silent = T,
      {
        subset <- hsc[,hsc$patient == i & (hsc$integrated_snn_res.0.05 == 1 | hsc$integrated_snn_res.0.05 == 2)]
        Idents(subset) <- "timing"
        
        HSCmarkers <- FindMarkers(subset, ident.1 = "post",
                                  group.by = "timing",
                                  verbose = T) %>%
          arrange(desc(avg_log2FC))
        
        write.table(HSCmarkers, col.names = NA, sep = "\t", file = paste0("../results/tables/deg-cluster12-", i, "-hs(p)c-post-vs-pre.tsv"))
        
        subset <- ScaleData(object = subset, features = rownames(HSCmarkers))
        
        topGenes <- rownames(rbind(head(HSCmarkers, 20), tail(HSCmarkers, 20)))

        ggsave(filename = paste0("../results/graphs/heatmap-cluster12-", i, "-hs(p)c-post-vs-pre.jpeg"),
               dpi = 300, height = 10, width = 8,
               plot = DoHeatmap(subset, group.by = "timing",
                                features = topGenes,
                                size = 2, raster = TRUE, label = FALSE) +
                 theme(axis.text.y = element_text(size = 7)))
      }
  )
  
}
#Trying to compare pre and post in only cluster 0 (nothing significant):
for (i in levels(aml$patient)) {
  try(silent = T,
      {
        subset <- hsc[,hsc$patient == i & hsc$integrated_snn_res.0.05 == 0]
        Idents(subset) <- "timing"
        
        HSCmarkers <- FindMarkers(subset, ident.1 = "post",
                                  group.by = "timing",
                                  verbose = T) %>% filter(p_val_adj < 0.05) %>%
          arrange(desc(avg_log2FC))
        
        write.table(HSCmarkers, col.names = NA, sep = "\t", file = paste0("../results/tables/deg-cluster0-", i, "-hs(p)c-post-vs-pre.tsv"))
        
        subset <- ScaleData(object = subset, features = rownames(HSCmarkers))
        
        topGenes <- rownames(rbind(head(HSCmarkers, 20), tail(HSCmarkers, 20)))

        ggsave(filename = paste0("../results/graphs/heatmap-cluster0-", i, "-hs(p)c-post-vs-pre.jpeg"),
               dpi = 300, height = 10, width = 8,
               plot = DoHeatmap(subset, group.by = "timing",
                                features = topGenes,
                                size = 2, raster = TRUE, label = FALSE) +
                 theme(axis.text.y = element_text(size = 7)))
      }
  )
  
}
