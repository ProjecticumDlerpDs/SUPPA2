#script filtering met seurat

#inladen packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)

#inladen data
features_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_feature_metadata.csv.gz")
samples_metadata <- read.csv("/Users/jalisavanderzeeuw/Desktop/DesktopJalisa/Rstudio/suppa2/seurat/data/e85_sample_metadata.csv")

