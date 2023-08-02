#This is an R script to create multiple Heatmaps on CLL Seurat objects, post-scFEA analysis focusing on Supermodules

library(Seurat)
library(readr)
library(tidyseurat)
library(ggplot2)

#Load Supermodules and Seurat object for post-scFEA analysis
CLL_d0_scFEA <- readRDS("E:/CLL rendeiro/CLL_d0_scFEA.rds")
CLL_d30_scFEA <- readRDS("E:/CLL rendeiro/CLL_d30_scFEA.rds")
CLL_d120_scFEA <- readRDS("E:/CLL rendeiro/CLL_d120_scFEA.rds")
CLL_d150_scFEA <- readRDS("E:/CLL rendeiro/CLL_d150_scFEA.rds")
CLL_d280_scFEA <- readRDS("E:/CLL rendeiro/CLL_d280_scFEA.rds")

Supermodules <- read_csv("E:/CLL rendeiro/Human_M168_information.symbols.csv")
# Supermodules
Glycolysis <- with(Supermodules, Modules[Supermodule_id == 1])
PPP <- with(Supermodules, Modules[Supermodule_id == 3])
FA <- with(Supermodules, Modules[Supermodule_id == 4])
Aspartate <- with(Supermodules, Modules[Supermodule_id == 5])
Glutamate <- with(Supermodules, Modules[Supermodule_id == 8])
Urea <- with(Supermodules, Modules[Supermodule_id == 10])
Transport <- with(Supermodules, Modules[Supermodule_id == 12])
Glycogen <- with(Supermodules, Modules[Supermodule_id == 14])
Purines <- with(Supermodules, Modules[Supermodule_id == 20])
Pyrimidines <- with(Supermodules, Modules[Supermodule_id == 21])
Serine <- with(Supermodules, Modules[Supermodule_id == 2])
BAla <- with(Supermodules, Modules[Supermodule_id == 6])
Steroid <- with(Supermodules, Modules[Supermodule_id == 22])
NLink <- with(Supermodules, Modules[Supermodule_id == 15])
OLink <- with(Supermodules, Modules[Supermodule_id == 16])
Sialic <- with(Supermodules, Modules[Supermodule_id == 17])
Chondro <- with(Supermodules, Modules[Supermodule_id == 18])
Heparan <- with(Supermodules, Modules[Supermodule_id == 19])
Hyaluronic <- with(Supermodules, Modules[Supermodule_id == 13])
Spermine <- with(Supermodules, Modules[Supermodule_id == 11])
Leucine <- with(Supermodules, Modules[Supermodule_id == 9])
Propanoyl <- with(Supermodules, Modules[Supermodule_id == 7])

# Lets limit the analysis to specific cell clusters:
CLL_d0_plot <- CLL_d0_scFEA %>% filter(labels %in% c("B cells", "CD4+ T cells", "CD8+ T cells",
                                                     "Monocytes", "NK cells"))

CLL_d30_plot <- CLL_d30_scFEA %>% filter(labels %in% c("B cells", "CD4+ T cells", "CD8+ T cells",
                                                       "Monocytes", "NK cells")) 

CLL_d120_plot <- CLL_d120_scFEA %>% filter(labels %in% c("B cells", "CD4+ T cells", "CD8+ T cells",
                                                         "Monocytes", "NK cells")) 

CLL_d150_plot <- CLL_d150_scFEA %>% filter(labels %in% c("B cells", "CD4+ T cells", "CD8+ T cells",
                                                         "Monocytes", "NK cells")) 

CLL_d280_plot <- CLL_d280_scFEA %>% filter(labels %in% c("B cells", "CD4+ T cells", "CD8+ T cells",
                                                         "Monocytes", "NK cells")) 

seurat_list <- c("D0_Pre_Ibr" = CLL_d0_plot, 
                 "D30_Post_Ibr"= CLL_d30_plot, 
                "CLL120_Post_Ibr" = CLL_d120_plot,
                "CLL_d150_Post_Ibr" = CLL_d150_plot, 
                "CLL_d280_Post_Ibr" = CLL_d280_plot)

feature_list <- list(
  "Glycolysis" = Glycolysis,
  "PPP" = PPP,
  "FA" = FA,
  "Aspartate" = Aspartate,
  "Glutamate" = Glutamate,
  "Urea" = Urea,
  "Transport" = Transport,
  "Glycogen" = Glycogen,
  "Purines" = Purines,
  "Pyrimidines" = Pyrimidines,
  "Serine" = Serine,
  "BAla" = BAla,
  "Steroid" = Steroid,
  "NLink" = NLink,
  "OLink" = OLink,
  "Sialic" = Sialic,
  "Chondro" = Chondro,
  "Heparan" = Heparan,
  "Hyaluronic" = Hyaluronic,
  "Spermine" = Spermine,
  "Leucine" = Leucine,
  "Propanoyl" = Propanoyl
)

generate_heatmap_pdf <- function(feature_list, output_file_prefix, data_frame_list, pdf_width = 8, pdf_height = 6) { 
  for (i in seq_along(data_frame_list)) {
    pdf(file = paste0(output_file_prefix, "_", names(data_frame_list)[i], ".pdf"), width = pdf_width, height = pdf_height)
    
    for (j in seq_along(feature_list)) {
      feature <- feature_list[[j]]
      feature_title <- names(feature_list)[j]  # Use the title from the feature_list
      
      feature_heatmap <- DoHeatmap(data_frame_list[[i]], features = feature, assay = 'FLUX', slot = 'scale.data', 
                                   group.by = 'labels', size = 4) + scale_fill_gradientn(colors = c("dodgerblue", "snow", "firebrick1")) + NoLegend() + ggtitle(paste(names(data_frame_list)[i], "-", feature_title))
      
      # You can customize the plot, such as changing the resolution, width, and height
      print(feature_heatmap)
    }
    
    dev.off()
  }
}



generate_heatmap_pdf(feature_list = feature_list, 
                       output_file_prefix = "heatmap_output", 
                       data_frame_list = seurat_list, 
                       pdf_width = 10, pdf_height = 8)



