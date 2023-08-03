#This is an R script to retrieve characteristic cluster modules from the CLL Seurat objects containing
#the "FLUX" assay (post-scFEA analysis)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyseurat)
library(patchwork)

#1. Load the Seurat objects
#Load Supermodules and Seurat object for post-scFEA analysis
CLL_d0_scFEA <- readRDS("E:/CLL rendeiro/CLL_d0_scFEA.rds")
CLL_d30_scFEA <- readRDS("E:/CLL rendeiro/CLL_d30_scFEA.rds")
CLL_d120_scFEA <- readRDS("E:/CLL rendeiro/CLL_d120_scFEA.rds")
CLL_d150_scFEA <- readRDS("E:/CLL rendeiro/CLL_d150_scFEA.rds")
CLL_d280_scFEA <- readRDS("E:/CLL rendeiro/CLL_d280_scFEA.rds")

#2. Calculate differentially expressed modules (one fluxomic cluster vs rest):
seurat_list <-  list("D0_scFEA" = CLL_d0_scFEA,"D30_scFEA" = CLL_d30_scFEA, "D120_scFEA" =CLL_d120_scFEA, 
                     "D150_scFEA" = CLL_d150_scFEA,"D280_scFEA" =CLL_d280_scFEA)

marker_dataframes <- list()

for (i in seq_along(seurat_list)) {
  markers <- FindAllMarkers(seurat_list[[i]], only.pos = FALSE, logfc.threshold = 0.1, 
                            features = VariableFeatures(object = seurat_list[[i]]), assay = 'FLUX', 
                            slot = "scale.data", test.use = "wilcox")
  custom_name <- paste(names(seurat_list)[i], "_Markers", sep = "")
  marker_dataframes[[custom_name]] <- as.data.frame(markers)
  filename <- paste0("Differentially_expressed_modules_", names(seurat_list)[i], ".csv")
  write.csv(markers, filename, row.names = TRUE)
}

top10_markers_list <- list()

# Perform the operations to get the top-10 markers for each cluster after the initial loop
for (i in seq_along(marker_dataframes)) {
  markers <- marker_dataframes[[i]]
  
  
  # Group by cluster and get the top-10 markers for each cluster
  top10_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_diff)
  
  # Store the top-10 marker dataframes for each cluster in the list
  top10_markers_list[[names(marker_dataframes)[i]]] <- as.data.frame(top10_markers)
}


p1 <- DoHeatmap(CLL_d0_scFEA, features = top10_markers_list$D0_scFEA_Markers$gene, assay = 'FLUX', slot = 'scale.data', 
                                    size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle("Modules for metabolic clusters, CLL_d0")

p2 <- DoHeatmap(CLL_d30_scFEA, features = top10_markers_list$D30_scFEA_Markers$gene, assay = 'FLUX', slot = 'scale.data', 
                size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle("Modules for metabolic clusters, CLL_d30")

p3 <- DoHeatmap(CLL_d120_scFEA, features = top10_markers_list$D120_scFEA_Markers$gene, assay = 'FLUX', slot = 'scale.data', 
                size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle("Modules for metabolic clusters, CLL_d120")

p4 <- DoHeatmap(CLL_d150_scFEA, features = top10_markers_list$D150_scFEA_Markers$gene, assay = 'FLUX', slot = 'scale.data', 
                size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle("Modules for metabolic clusters, CLL_d150")

p5 <- DoHeatmap(CLL_d280_scFEA, features = top10_markers_list$D280_scFEA_Markers$gene, assay = 'FLUX', slot = 'scale.data', 
                      size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle("Modules for metabolic clusters, CLL_d280")
