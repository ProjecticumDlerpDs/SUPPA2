#script filtering met seurat

#inladen packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(Matrix)

#inladen data
features_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_feature_metadata.csv.gz")
samples_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_sample_metadata.csv")

#inladen MTX file
counts <- ReadMtx(
  mtx = "/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/e85_count_matrix.mtx.gz",                 
  features = "/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_feature_metadata.csv.gz",          
  cells =       "/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_sample_metadata.csv",
  feature.sep = ",",
  feature.column = 1,
  cell.sep = ",",
  cell.column = 3,
  skip.cell = 1,
  skip.feature = 1
)

#aanmaken seurat object met ruwe data
seurat <- CreateSeuratObject(counts = counts,
                             project = "mouse_embryo",
                             min.cells = 3,
                             min.features = 200)

#bekijk samenvatting en check of de eerste paar rijen kloppen
seurat
head(seurat)

#bereken percent.mt
mito.features <- grep("^mt-", rownames(seurat), value = TRUE, ignore.case = TRUE)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mito.features)

#QC (features, countes, percent.mt)
png("seurat/output/qc-vlnplot.png", width = 800, height = 600)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        layer = "counts",
        ncol = 3)

dev.off()

#QC scatternplot 
png("seurat/output/scatter_plots.png", width = 1200, height = 800)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
print(plot1 + plot2)

dev.off()

#QC filtering, subset van de data
seurat.filtered <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 6)

#normalisatie, selectie variabele genen, regressie op mitochondriale genexpressie & dataschaling
seurat.filtered <- SCTransform(seurat.filtered, vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(seurat.filtered) <- "SCT"

#toon hoeveel variabele genen er geselecteerd zijn
length(VariableFeatures(seurat.filtered))

#10 meest variabele genen uitzoeken 
top10 <- head(VariableFeatures(seurat.filtered), 10)
print(top10)

#visualisatie van de hvg's
png("seurat/output/variable_features_plot.png", width = 1000, height = 600)
plot1 <- VariableFeaturePlot(seurat.filtered)
top10 <- head(VariableFeatures(seurat.filtered), 10)
print(plot1 + LabelPoints(plot1, points = top10, repel = TRUE))

dev.off()

#PCA analyse
seurat.filtered <- RunPCA(seurat.filtered, features = VariableFeatures(object = seurat.filtered))
print(seurat.filtered[["pca"]], dims = 1:5, nfeatures = 5)

#sla loadings plot op als afbeelding
png("seurat/output/PCA_loadings_plot.png", width = 800, height = 600)
VizDimLoadings(seurat.filtered, dims = 1:2, reduction = "pca")

dev.off()
#Deze plots helpen de belangrijkste genen te identificeren die de variatie in de dataset verklaren.
#Dit is nuttig om de biologie achter de data te begrijpen, de kwaliteit van de data te checken,
# en context te geven aan de splicing-analyse met suppa2.

#visualisatie van cellen in PC-ruimte
png("seurat/output/PCA_plot.png", width = 800, height = 600)
DimPlot(seurat.filtered, reduction = "pca")

dev.off()

#PCA analyse plotten 
png("seurat/output/elbow_plot.png", width = 800, height = 600)
ElbowPlot(seurat.filtered)
dev.off()

#vind buren op basis van PC1 t/m PC8 en cluster de cellen
seurat.filtered <- FindNeighbors(seurat.filtered, dims = 1:8)
seurat.filtered <- FindClusters(seurat.filtered,resolution = 0.5)

#visualiseer clusters met UMAP
seurat.filtered <- RunUMAP(seurat.filtered, dims = 1:8)
umap_plot <- DimPlot(seurat.filtered, reduction = "umap", label = TRUE)
print(umap_plot)

#opslaan van de UMAP plot
png("seurat/output/umap_clusters.png", width = 800, height = 600)
print(umap_plot)

dev.off()


#bekijk top 3 markergenen per cluster
markers_all <- FindAllMarkers(seurat.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top3_markers_per_cluster <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3)
print(top3_markers_per_cluster, n = Inf)

# Gekozen voor cluster 2 (splicing-actief) en cluster 3 (stabiel en homogeen)

#vergelijk clusters 2 en 3
top_genen <- markers_all %>%
  filter(cluster %in% c(2, 3)) %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

#top15 per cluster apart
top_2 <- markers_all %>% filter(cluster == 2) %>% top_n(15, avg_log2FC) %>% pull(gene)
top_3 <- markers_all %>% filter(cluster == 3) %>% top_n(15, avg_log2FC) %>% pull(gene)

# Overlap bekijken tussen cluster 2 en 3
overlap <- intersect(top_2, top_3)
print(overlap)

#maak subset van alleen cluster 2 en 3
seurat_2_3 <- subset(seurat.filtered, idents = c(2,3))

#combineer de top genen van beide clusters
top_genes <- c(top_2, top_3) %>% unique()

#maak dotplot voor alleen die clusters met die genen
DotPlot(seurat_2_3, features = top_genes) + RotatedAxis()

png("seurat/output/dotplot_2_3.png", width = 800, height = 600)
print(DotPlot(seurat_2_3, features = top_genes) + RotatedAxis())
dev.off()

