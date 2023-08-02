# This is the R script to prepare integrate scFEA outputs regarding cell-wise fluxomics in the CLL Seurat objects created before:

#1. Load scFEA outputs
#Day 0
flux0 <- read.csv('CLL_grouped_D0_flux.csv', header = T, row.names = 1)
flux0 <- data.matrix(flux0)
flux0 <- t(flux0)
#Day 30
flux30 <- read.csv('CLL_grouped_D30_flux.csv', header = T, row.names = 1)
flux30 <- data.matrix(flux30)
flux30 <- t(flux30)
#Day 120
flux120 <- read.csv('CLL_grouped_D120_flux.csv', header = T, row.names = 1)
flux120 <- data.matrix(flux120)
flux120 <- t(flux120)
#Day 150
flux150 <- read.csv('CLL_grouped_D150_flux.csv', header = T, row.names = 1)
flux150 <- data.matrix(flux150)
flux150 <- t(flux150)
#Day 280
flux280 <- read.csv('CLL_grouped_D280_flux.csv', header = T, row.names = 1)
flux280 <- data.matrix(flux280)
flux280 <- t(flux280)


#2. Integrate flux files as a new assay in the respective Seurat objects:
CLL0[["FLUX"]] <- CreateAssayObject(counts = flux0)
DefaultAssay(CLL0) <- 'FLUX'

CLL30[["FLUX"]] <- CreateAssayObject(counts = flux30)
DefaultAssay(CLL30) <- 'FLUX'

CLL120[["FLUX"]] <- CreateAssayObject(counts = flux120)
DefaultAssay(CLL120) <- 'FLUX'

CLL150[["FLUX"]] <- CreateAssayObject(counts = flux150)
DefaultAssay(CLL150) <- 'FLUX'

CLL280[["FLUX"]] <- CreateAssayObject(counts = flux280)
DefaultAssay(CLL280) <- 'FLUX'



# 3. Fluxomic clustering of Seurat objects
Srt_list <- list(CLL0, CLL30, CLL120, CLL150, CLL280)
for (i in 1:length(Srt_list)) {
  Srt_list[[i]] <- FindVariableFeatures(Srt_list[[i]],  selection.method = "vst", nfeatures = 2000, verbose = T, assay = "FLUX")
  Srt_list[[i]] <- ScaleData(Srt_list[[i]], features = rownames(Srt_list[[i]]), assay = 'FLUX', verbose = TRUE)
  Srt_list[[i]] <- RunPCA(Srt_list[[i]], features = VariableFeatures(object = Srt_list[[i]]), npcs = 10, 
                reduction.name = 'pca.flux', verbose = T, assay = "FLUX")
  Srt_list[[i]] <- FindNeighbors(Srt_list[[i]], dims = 1:10, verbose = F, reduction = "pca.flux")
  Srt_list[[i]] <- FindClusters(Srt_list[[i]], resolution = 0.2, verbose = F)
  Srt_list[[i]] <- RunUMAP(Srt_list[[i]], dims = 1:10, assay = 'FLUX', reduction.name = "umap.flux", verbose = T, reduction = "pca.flux")
}

for (i in 1:length(Srt_list)) {                                          #sanity check
  print(DimPlot(Srt_list[[i]], reduction = "umap.flux", group.by = "first.labels") + ggtitle('UMAP of Flux with cell annotation'))
  print(DimPlot(Srt_list[[i]], reduction = "umap.flux", label = T) + ggtitle('UMAP of Flux with scFEA clusters'))
}


#Split list
CLL_d0_scFEA <- Srt_list[[1]]
CLL_d30_scFEA <- Srt_list[[2]]
CLL_d120_scFEA <- Srt_list[[3]]
CLL_d150_scFEA <- Srt_list[[4]]
CLL_d280_scFEA <- Srt_list[[5]]

saveRDS(CLL_d0_scFEA, "CLL_d0_scFEA.rds")
saveRDS(CLL_d30_scFEA, "CLL_d30_scFEA.rds")
saveRDS(CLL_d120_scFEA, "CLL_d120_scFEA.rds")
saveRDS(CLL_d150_scFEA, "CLL_d150_scFEA.rds")
saveRDS(CLL_d280_scFEA, "CLL_d280_scFEA.rds")
