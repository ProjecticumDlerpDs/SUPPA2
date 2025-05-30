#script filtering met seurat

#inladen packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)

#inladen data
features_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/Desktop - Jalisa/Rstudio/suppa2/seurat/data/feature_metadata.csv")
samples_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/Desktop - Jalisa/Rstudio/suppa2/seurat/data/e85_sample_metadata.csv")

