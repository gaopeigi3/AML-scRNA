#################Environment####################################################
#Loading libraries:
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(patchwork)
library(GenomicRanges)

#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Enable parallelization
plan("multisession", workers = 24)

#Clearing memory:
gc()

aml <- readRDS("../data/byproducts/01-aml-combined.rds")

#################Quality Control################################################

##Annotating peaks:
#Getting EnsDb annotations:
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = T)

#Changing annotation style to UCSC:
seqlevelsStyle(annotations) <- "UCSC"
Signac::genome(annotations) <- "hg38"

#Adding annotations to seurat object:
Annotation(aml) <- annotations

##Nucleosome and TSS enrichment:
#Nucleosome signal score/cell:
aml <- NucleosomeSignal(aml)

#TSS enrichment score/cell:
aml <- TSSEnrichment(aml, fast = F)

#Adding blacklist ration and reads/peak:
aml$pct_reads_in_peaks <- aml$peak_region_fragments / aml$passed_filters * 100
aml$blacklist_ratio <- aml$blacklist_region_fragments / aml$peak_region_fragments

#Plotting TSS enrichment:
aml$high.tss <- ifelse(aml$TSS.enrichment > 2, "High", "Low")
ggsave(filename = "../results/graphs/tss-enrichment-by-level.jpg",
       width = 10, height = 7, dpi = 300,
       plot = TSSPlot(aml, group.by = "high.tss") + NoLegend())
ggsave(filename = "../results/graphs/tss-enrichment-by-sample.jpg",
       width = 10, height = 7, dpi = 300,
       plot = TSSPlot(aml, group.by = "sample") + NoLegend())

#Plotting nucleosome banding:
aml$nucleosome_group <- ifelse(aml$nucleosome_signal > 4, "NS > 4", "NS < 4")
ggsave(filename = "../results/graphs/nucleosome-banding.jpg",
       width = 10, height = 7, dpi = 300,
       plot = FragmentHistogram(aml, group.by = "nucleosome_group"))


#Adding total peak fragment distribution:
aml$peak_region_fragments_log <- log2(aml$peak_region_fragments + 1)
#Plotting the above:
p1 <- aml@meta.data %>%
    ggplot(aes(x = peak_region_fragments)) +
    geom_histogram(bins = 100, col = "black", fill = "grey60") +
    theme_classic() +
    geom_vline(xintercept = c(500, 20000), col = "red") +
    ggtitle("Total peak fragments histogram")

p2 <- aml@meta.data %>%
    ggplot(aes(x = peak_region_fragments_log)) +
    geom_histogram(bins = 100, col = "black", fill = "grey60") +
    theme_classic() +
    geom_vline(xintercept = c(log2(500 + 1), log2(20000+1)), col = "red") +
    ggtitle("Log2 total peak fragments histogram")

ggsave(filename = "../results/graphs/peak-fragment.jpg",
       width = 10, height = 7, dpi = 300,
       plot = p1 + p2)

#Scatter plots:
p.list <- lapply(unique(aml$sample), function(sample){
  FeatureScatter(aml[, aml$sample == sample], "peak_region_fragments_log", "TSS.enrichment") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_vline(xintercept = 3, linetype = "dashed") +
    expand_limits(x = 0, y = 0) +
    ggtitle(sample) +
    NoLegend()
})
ggsave(filename = "../results/graphs/fragment-tss-enrichment.jpg",
       width = 15, height = 15, dpi = 300,
       plot = wrap_plots(p.list))

VlnPlot(
  object = aml,
  features = c('pct_reads_in_peaks', 'peak_region_fragments_log',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
) *
    geom_boxplot(width = .3, fill = "white")


#Filtering the object:
aml <- subset(aml, subset = 
                peak_region_fragments > 500 &
                peak_region_fragments < 20000 &
                pct_reads_in_peaks > 25 &
                blacklist_ratio < 0.05 &
                nucleosome_signal < 4 &
                TSS.enrichment > 2)

#Violin plots after filtration:
ggsave(filename = "../results/graphs/general-qc-metrics.jpg",
       width = 15, height = 10, dpi = 300,
       plot = VlnPlot(
         object = aml,
         features = c('pct_reads_in_peaks', 'peak_region_fragments_log',
                      'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
         pt.size = 0,
         ncol = 5
       ) *
         geom_boxplot(width = .3, fill = "white"))

saveRDS(aml, "../data/byproducts/02-aml-filtered.rds", compress = F)

#################Add Metadata###################################################

#Creating sample column:
aml$sample <- sub("(.*?)_{1}(.*?)($|-.*)", "\\1", colnames(aml))
aml$sample <- factor(aml$sample, levels = c("Healthy1_control", "Healthy2_control",
                                            "Patient1_pre", "Patient1_post", "Patient1_post2",
                                            "Patient2_pre", "Patient2_post"))

#Creating patient column:
aml$patient <- sub("(.*)_(.*)", "\\1", aml$sample)
aml$patient <- sub("([[:alpha:]])([[:digit:]])", "\\1 \\2", aml$patient)
aml$patient <- factor(aml$patient, levels = c("Healthy 1", "Healthy 2",
                                              "Patient 1", "Patient 2"))

#Creating timing column:
aml$timing <- sub("(.*)_(.*)", "\\2", aml$sample)
aml$timing <- factor(aml$timing, levels = c("pre", "post"))

#Creating collection date column:
aml$collectionDate <- NA_Date_
aml$collectionDate[aml$sample == "Patient1_pre"] <- "2018-05-17"
aml$collectionDate[aml$sample == "Patient1_post"] <- "2018-07-17"
aml$collectionDate[aml$sample == "Patient1_post2"] <- "2018-08-30"
aml$collectionDate[aml$sample == "Patient2_pre"] <- "2018-07-27"
aml$collectionDate[aml$sample == "Patient2_post"] <- "2018-09-05"
aml$collectionDate <- as.Date(aml$collectionDate)

#Creating blast column:
aml$blast <- NA_integer_
aml$blast[aml$sample == "Patient1_pre"] <- 0.34
aml$blast[aml$sample == "Patient1_post"] <- 0.09
aml$blast[aml$sample == "Patient1_post2"] <- 0.50
aml$blast[aml$sample == "Patient2_pre"] <- 0.67
aml$blast[aml$sample == "Patient2_post"] <- 0.40

#Creating category column:
aml$category <- "healthy"
aml$category[aml$patient == "Patient 1"] <- "R/R"
aml$category[aml$patient == "Patient 2"] <- "R/R"

#Creating response column:
aml$response <- "healthy"
aml$response[aml$patient == "Patient 1"] <- "Relapse"
aml$response[aml$patient == "Patient 2"] <- "Refractory"

saveRDS(aml, "../data/byproducts/02-aml-filtered.rds", compress = F)