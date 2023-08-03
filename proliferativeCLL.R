DimPlot(CLL_d0_scFEA, label = T)

VlnPlot(CLL_d0_scFEA, features = c("CXCR4", "CD5"), assay = "SCT")
FeaturePlot(CLL_d0_scFEA, features = c("CXCR4", "CD5"), assay = "SCT", blend = T)


Nebulosa::plot_density(CLL_d0_scFEA, features = c("CXCR4", "CD5"), assay="SCT",
                       joint = TRUE, pal = "viridis")



Idents(CLL_d0_scFEA)
CLL_d0_scFEA$Scfea_clustering <- Idents(CLL_d0_scFEA)



DefaultAssay(CLL_d0_scFEA) <- "SCT"

DimPlot(CLL_d0_scFEA)
FeaturePlot(CLL_d0_scFEA, features = c("CXCR4", "CD5"), blend = T)
Nebulosa::plot_density(CLL_d0_scFEA, features = c("CXCR4", "CD5"),
                       joint = TRUE, pal = "viridis", reduction = "umap")
FeaturePlot(CLL_d0_scFEA, features = c("CXCR4", "CD5"))

FeaturePlot(CLL_d0_scFEA, features = 'CD5')

# (CLL_d0_scFEA@assays$SCT@counts["CD5",]

# Replace 'gene_of_interest' with the actual name of the gene you want to filter for
gene_of_interest <- "CD5"

# Set a threshold for the expression level to define positive cells
DefaultAssay(CLL_d0_scFEA) <- 'SCT'
seurat_filtered <- Seurat::WhichCells(CLL_d0_scFEA, expression = CD5 > 0.1 & CXCR4 > 0.1)

filter <- CLL_d0_scFEA[, seurat_filtered]
filter <- filter %>% tidyseurat::filter(labels %in% c('B cells', 'Basophils'))

DimPlot(filter, group.by = "labels")
DimPlot(CLL_d0_scFEA)

FeaturePlot(filter, features = c('CD5', 'CXCR4'), blend = T)
Nebulosa::plot_density(filter, features = c("CD5", "CXCR4"),
                       joint = TRUE, pal = "viridis", reduction = 'umap')

DimPlot(filter, reduction = 'umap.flux')
FeaturePlot(filter, features = c('CD5', 'CXCR4', 'CD38'), reduction = 'umap.flux')
Nebulosa::plot_density(filter, features = c("CD5", "CD38"),
                       joint = TRUE, pal = "viridis", reduction = 'umap.flux')

Nebulosa::plot_density(filter, features = c("CXCR4"),
                       joint = TRUE, pal = "viridis", reduction = 'umap.flux')


p1 <- Nebulosa::plot_density(filter, features = c("CD5"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p2 <- Nebulosa::plot_density(filter, features = c("CD38"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p3 <- Nebulosa::plot_density(filter, features = c("CXCR4"),
                             joint = TRUE, pal = "viridis", reduction = 'umap.flux')

p4 <- DimPlot(filter, reduction = 'umap.flux', label = T, group.by = 'Scfea_clustering')

p1 | p2 | p3 | p4


p5 <- Nebulosa::plot_density(filter, features = c("NFKB1"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p6 <- Nebulosa::plot_density(filter, features = c("NOTCH1"),
                             joint = FALSE, pal = "viridis", reduction = 'umap.flux')
p7 <- Nebulosa::plot_density(filter, features = c("NFATC1"),
                             joint = TRUE, pal = "viridis", reduction = 'umap.flux')
p8 <- Nebulosa::plot_density(filter, features = c("FOXC1"),
                             joint = TRUE, pal = "viridis", reduction = 'umap.flux')

p5 | p6 | p7 




#filter <- CellCycleScoring(filter, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
filter <- SCTransform(filter, vst.flavor = "v2", verbose = FALSE)
filter <- RunPCA(filter, features=VariableFeatures(filter))

# filter <- RunUMAP(filter, dims=1:30, reduction = "pca", verbose = FALSE, n.epochs = 200, min.dist = 0.1) #, n.neighbors = 20)
nPCs = 30
res = 0.5
filter <- FindNeighbors(filter, dims=1:nPCs)
filter <- FindClusters(filter, resolution=res)

filter <- RunUMAP(filter, dims=1:30, reduction = "pca", verbose = FALSE) #, n.epochs = 200, min.dist = 0.1, n.neighbors = 20)

p9 <- DimPlot(filter,  label = TRUE)
p10 <- DimPlot(filter, group.by = "labels", label = TRUE)
p11 <- DimPlot(filter, reduction = 'umap.flux')
p9 + p10 + p11 + plot_annotation(tag_levels = 'A')



p1 <- Nebulosa::plot_density(filter, features = c("CD5"),
                             joint = FALSE, pal = "viridis", reduction = 'umap')
p2 <- Nebulosa::plot_density(filter, features = c("CD38"),
                             joint = FALSE, pal = "viridis", reduction = 'umap')
p3 <- Nebulosa::plot_density(filter, features = c("CXCR4"),
                             joint = TRUE, pal = "viridis", reduction = 'umap')

p4 <- DimPlot(filter, reduction = 'umap.flux', label = T)

p1 | p2 | p3 | p5 | p7



p5 <- Nebulosa::plot_density(filter, features = c("NFKB1"),
                             joint = FALSE, pal = "viridis", reduction = 'umap')
p6 <- Nebulosa::plot_density(filter, features = c("NOTCH1"),
                             joint = FALSE, pal = "viridis", reduction = 'umap')
p7 <- Nebulosa::plot_density(filter, features = c("NFATC1"),
                             joint = TRUE, pal = "viridis", reduction = 'umap')
p8 <- Nebulosa::plot_density(filter, features = c("FOXC1"),
                             joint = TRUE, pal = "viridis", reduction = 'umap')

p5 | p7 

p12 <- Nebulosa::plot_density(filter, features = c("CD5", "CD38", "NFATC1"),
                              joint = TRUE, pal = "viridis", reduction = 'umap')
p12

Idents(filter)
f.markers <- FindAllMarkers(filter, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = 0.25)
f.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

p10 <- DoHeatmap(filter, features = top10$gene) 
p10

p11 <- DoHeatmap(filter, features = c("CD5", "CD38", "NFATC1"), slot = 'counts')
p11


write.csv(f.markers, "f_markers.csv")
