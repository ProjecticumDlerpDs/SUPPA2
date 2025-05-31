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
png("qc-vlnplot.png", width = 800, height = 600)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        layer = "counts",
        ncol = 3)

dev.off()

#QC scatternplot 
png("scatter_plots.png", width = 1200, height = 800)
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
png("variable_features_plot.png", width = 1000, height = 600)
plot1 <- VariableFeaturePlot(seurat.filtered)
top10 <- head(VariableFeatures(seurat.filtered), 10)
print(plot1 + LabelPoints(plot1, points = top10, repel = TRUE))

dev.off()

#PCA analyse
seurat.filtered <- RunPCA(seurat.filtered, features = VariableFeatures(object = seurat.filtered))
print(seurat.filtered[["pca"]], dims = 1:5, nfeatures = 10)

#sla loadings plot op als afbeelding
png("PCA_loadings_plot.png", width = 800, height = 600)
VizDimLoadings(seurat.filtered, dims = 1:2, reduction = "pca")

dev.off()
#Deze plots helpen de belangrijkste genen te identificeren die de variatie in de dataset verklaren.
#Dit is nuttig om de biologie achter de data te begrijpen, de kwaliteit van de data te checken,
# en context te geven aan de splicing-analyse met suppa2.

#visualisatie van cellen in PC-ruimte
png("PCA_plot.png", width = 800, height = 600)
DimPlot(seurat.filtered, reduction = "pca")

dev.off()
