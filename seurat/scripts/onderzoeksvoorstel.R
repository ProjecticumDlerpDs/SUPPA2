# umap die staat in Seurat.Rmd
seurat.filtered <- FindNeighbors(seurat.filtered, dims = 1:8)
seurat.filtered <- FindClusters(seurat.filtered, resolution = 0.5)
seurat.filtered <- RunUMAP(seurat.filtered, dims = 1:8)
DimPlot(seurat.filtered, reduction = "umap", label = TRUE)

# nieuwe clustering met een andere resolutie, in een nieuw object
seurat.new <- seurat.filtered

# Nieuwe clustering met resolutie 0.5 bijvoorbeeld
seurat.new <- FindNeighbors(seurat.new, dims = 1:8)
seurat.new <- FindClusters(seurat.new, resolution = 1.5)
seurat.new <- RunUMAP(seurat.new, dims = 1:8)

# Visualiseren
DimPlot(seurat.new, reduction = "umap", label = TRUE)

# rint aantal cellen per cluster
cluster_sizes <- table(Idents(seurat.filtered))
print(cluster_sizes)



