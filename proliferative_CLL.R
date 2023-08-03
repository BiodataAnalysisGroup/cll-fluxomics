#This is an R script to locate where the metabolically active cluster of CLL clones in the peripheral blood resides in terms of 
# scRNA-seq metabolic clusters
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyseurat)
library(patchwork)
library (Nebulosa)
library(cowplot)

#1. Load the Seurat objects
#Load Supermodules and Seurat object for post-scFEA analysis
CLL_d0_scFEA <- readRDS("E:/CLL rendeiro/CLL_d0_scFEA.rds")
CLL_d30_scFEA <- readRDS("E:/CLL rendeiro/CLL_d30_scFEA.rds")
CLL_d120_scFEA <- readRDS("E:/CLL rendeiro/CLL_d120_scFEA.rds")
CLL_d150_scFEA <- readRDS("E:/CLL rendeiro/CLL_d150_scFEA.rds")
CLL_d280_scFEA <- readRDS("E:/CLL rendeiro/CLL_d280_scFEA.rds")

#2. Use Nebulosa to locate CLL cells that co-express CD5 and CD38, CXCR4, NFkB, NFACT1, NOTCH1
DefaultAssay(CLL_d0_scFEA) <- "SCT"
DefaultAssay(CLL_d30_scFEA) <- "SCT"
DefaultAssay(CLL_d120_scFEA) <- "SCT"
DefaultAssay(CLL_d150_scFEA) <- "SCT"
DefaultAssay(CLL_d280_scFEA) <- "SCT"

p1 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CD5"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p2 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CD38"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p3 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CXCR4"),
                             joint = TRUE, pal = "viridis", reduction = 'umap.flux')
p1 | p2 | p3


p3 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CD5"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p4 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CD38"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p5 <- Nebulosa::plot_density(CLL_d0_scFEA, features = c("CXCR4"),
                             joint = TRUE, pal = "viridis", reduction = 'umap.flux')
p3 | p4 | p5



seurat_list <- list("CLL_d0" = CLL_d0_scFEA, 
                    "CLL_d30" = CLL_d30_scFEA, 
                    "CLL_d120"= CLL_d120_scFEA, 
                    "CLL_d150" = CLL_d150_scFEA, 
                    "CLL_280" = CLL_d280_scFEA)




# Function to create the UMAP plots and save them to a PDF file
create_combined_umap_plots <- function(seurat_object, filename, pdf_width = 25, pdf_height = 8) {
  p1 <- Nebulosa::plot_density(seurat_object, features = c("CD5"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p2 <- Nebulosa::plot_density(seurat_object, features = c("CD38"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p3 <- Nebulosa::plot_density(seurat_object, features = c("CXCR4"),
                               joint = TRUE, pal = "viridis", reduction = 'umap.flux')
  
  # Combine the three plots side by side using the cowplot library
  combined_plot <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
  
  # Save the combined plot to a PDF file
  pdf(filename)
  print(combined_plot)
  dev.off()
}

# Loop over the list of Seurat objects and create PDF files for each one
for (i in seq_along(seurat_list)) {
  filename <- paste0("umap_plots_", names(seurat_list)[i], ".pdf")
  create_combined_umap_plots(seurat_list[[i]], filename, pdf_width = 25, pdf_height = 8)
}




# Function to create the UMAP plots and save them to a PDF file
create_combined_umap_plots <- function(seurat_object, filename) {
  p1 <- Nebulosa::plot_density(seurat_object, features = c("CD5"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p2 <- Nebulosa::plot_density(seurat_object, features = c("CD38"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p3 <- Nebulosa::plot_density(seurat_object, features = c("CXCR4"),
                               joint = TRUE, pal = "viridis", reduction = 'umap.flux')
  
  # Combine the three plots side by side using the cowplot library
  combined_plot <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
  
  # Save the combined plot to a PDF file with increased width
  pdf(filename, width = 15, height = 8)  # Adjust the width as needed (e.g., width = 15)
  print(combined_plot)
  dev.off()
}

# Loop over the list of Seurat objects and create PDF files for each one
for (i in seq_along(seurat_list)) {
  filename <- paste0("umap_plots_", names(seurat_list)[i], ".pdf")
  create_combined_umap_plots(seurat_list[[i]], filename)
}




# Function to create the UMAP plots and save them to a PDF file
create_combined_umap_plots2 <- function(seurat_object, filename) {
  p1 <- Nebulosa::plot_density(seurat_object, features = c("NFKB1"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p2 <- Nebulosa::plot_density(seurat_object, features = c("NOTCH1"),
                               joint = FALSE, pal = "viridis", reduction = 'umap.flux')
  p3 <- Nebulosa::plot_density(seurat_object, features = c("NFATC1"),
                               joint = TRUE, pal = "viridis", reduction = 'umap.flux')
  
  # Combine the three plots side by side using the cowplot library
  combined_plot <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
  
  # Save the combined plot to a PDF file with increased width
  pdf(filename, width = 15, height = 8)  # Adjust the width as needed (e.g., width = 15)
  print(combined_plot)
  dev.off()
}

# Loop over the list of Seurat objects and create PDF files for each one
for (i in seq_along(seurat_list)) {
  filename <- paste0("umap_plots_proliferative_genes_", names(seurat_list)[i], ".pdf")
  create_combined_umap_plots2(seurat_list[[i]], filename)
}
