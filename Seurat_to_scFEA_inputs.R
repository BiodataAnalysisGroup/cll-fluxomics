# This is the R script to prepare the necessary data from the Seurat objects from the Rendeiro study
# to provide it as input into scFEA neural networks.

#0. Load necessary libraries
library(Seurat)
library(tidyseurat)
library(ggplot2)
library(patchwork)

#1. Load the integrated Seurat  CLL object created in the previous script and subset to each give timepoint:
CLL <- readRDS("E:/cll_integrated_annot.rds") #use the path for your working directory
#Explore the Seurat object
DimPlot(CLL, group.by = "labels", split.by = "days")
DimPlot(CLL, group.by = "labels", split.by = "subject")
# Subset to each Day
CLL0 <- CLL %>% filter(days == 'd0')
CLL30 <- CLL %>% filter(days == 'd30')
CLL120 <- CLL %>% filter(days == 'd120')
CLL150 <- CLL %>% filter(days == 'd150')
CLL280 <- CLL %>% filter(days == 'd280')
seurat_init_list <- c("Day0" = CLL0, 
                      "Day30" = CLL30, 
                      "Day120"=CLL120, 
                      "Day150"= CLL150, 
                      "Day280" = CLL280)

#2. Export the necessary data from Seurat objects to be used by scFEA as inputs:
# a. Expression matrix (we choose from the SCT assay)
# b. Information about Idents
#a.
for (i in seq_along(seurat_init_list)) {
  Data <- seurat_init_list[[i]]@assays$SCT@data
  # Generate a unique CSV filename for each Seurat object (e.g., Data_CLL_D0_1.csv, Data_CLL_D0_2.csv, ...)
  filename <- paste0("Expression Matrix of", i, ".csv")
  # Write the data to a CSV file
  write.csv(Data, filename, row.names = TRUE)
}  

#b.
# Assuming you have a list of Seurat objects called seurat_list

# Loop through each Seurat object in the list
for (i in seq_along(seurat_init_list)) {
  # Perform the operation on the current Seurat object
  idents <- Idents(seurat_init_list[[i]])
  
  # Create a unique filename for each Seurat object (e.g., Idents_CLL_d0_1.RData, Idents_CLL_d0_2.RData, ...)
  filename <- paste0("Idents", "_", names(seurat_init_list)[i], ".RData")
  
  # Save the idents object to the RData file
  save(list = "idents", file = filename)
}
