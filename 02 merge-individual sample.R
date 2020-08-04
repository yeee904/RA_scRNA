## 
### ---------------
###
### Create: Shanzhao Jin
### Date: 2019/10/01
### Update Log: 2020/08/04  
### This script aims to merge all individual Seurat objects from the same tissue (PBMC or synovial tissue)
###
### ---------------


#/User/bin/env Rscript
## Seurat standard workplow
library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)
library(umap)
library(RColorBrewer)
library(cowplot)
library(randomcoloR)
library(corrplot)
library(DoubletFinder)

cell_type_cols <- c(brewer.pal(9, "Set1"),
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4",
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030",
                    "#7CFC00","#000000","#708090")
HC1_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/HC1-PBMC/GRCh38")
HC2_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/HC2-PBMC/GRCh38")
HC3_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/HC3-PBMC/GRCh38")
HC4_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/HC4-PBMC/GRCh38")
RA1_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA1-PBMC/GRCh38")
RA2_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA2-PBMC/GRCh38")
RA3_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA3-PBMC/GRCh38")
RA4_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA04-PBMC/GRCh38/")
RA5_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA05-PBMC/GRCh38/")
RA6_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA06-PBMC/GRCh38/")
RA8_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA08-PBMC/GRCh38")
RA9_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA09-PBMC/GRCh38")
RA10_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA10-PBMC/GRCh38")
RA11_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA11-PBMC/GRCh38")
RA12_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA12-PBMC/GRCh38")
RA13_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA13-PBMC/GRCh38/")
RA14_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA14-PBMC/GRCh38/")
RA15_PBMC.data <- Read10X(data.dir = "D:/RA/samples/PBMC/RA15-PBMC/GRCh38")
# Singe Seurat object built -----------------------------------------------
# Set up HC1_PBMC_seurat object
HC1_PBMC_seurat <- CreateSeuratObject(counts = HC1_PBMC.data, project = "HC1_PBMC", min.cells = 5, min.features = 500)
HC1_PBMC_seurat$samples <- "HC1_PBMC"
HC1_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(HC1_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_HC1_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(HC1_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HC1_PBMC_seurat <- subset(HC1_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
HC1_PBMC_seurat <- NormalizeData(HC1_PBMC_seurat, verbose = FALSE)
HC1_PBMC_seurat <- FindVariableFeatures(HC1_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(HC1_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up HC2_PBMC_seurat object
HC2_PBMC_seurat <- CreateSeuratObject(counts = HC2_PBMC.data, project = "HC2_PBMC", min.cells = 5, min.features = 500)
HC2_PBMC_seurat$samples <- "HC2_PBMC"
HC2_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(HC2_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_HC2_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(HC2_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HC2_PBMC_seurat <- subset(HC2_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
HC2_PBMC_seurat <- NormalizeData(HC2_PBMC_seurat, verbose = FALSE)
HC2_PBMC_seurat <- FindVariableFeatures(HC2_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(HC2_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up HC3_PBMC_seurat object
HC3_PBMC_seurat <- CreateSeuratObject(counts = HC3_PBMC.data, project = "HC3_PBMC", min.cells = 5, min.features = 500)
HC3_PBMC_seurat$samples <- "HC3_PBMC"
HC3_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(HC3_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_HC3_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(HC3_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HC3_PBMC_seurat <- subset(HC3_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
HC3_PBMC_seurat <- NormalizeData(HC3_PBMC_seurat, verbose = FALSE)
HC3_PBMC_seurat <- FindVariableFeatures(HC3_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(HC3_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up HC4_PBMC_seurat object
HC4_PBMC_seurat <- CreateSeuratObject(counts = HC4_PBMC.data, project = "HC4_PBMC", min.cells = 5, min.features = 500)
HC4_PBMC_seurat$samples <- "HC4_PBMC"
HC4_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(HC4_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_HC4_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(HC4_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HC4_PBMC_seurat <- subset(HC4_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
HC4_PBMC_seurat <- NormalizeData(HC4_PBMC_seurat, verbose = FALSE)
HC4_PBMC_seurat <- FindVariableFeatures(HC4_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(HC4_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA1_PBMC_seurat object
RA1_PBMC_seurat <- CreateSeuratObject(counts = RA1_PBMC.data, project = "RA1_PBMC", min.cells = 5, min.features = 500)
RA1_PBMC_seurat$samples <- "RA1_PBMC"
RA1_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA1_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA1_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA1_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA1_PBMC_seurat <- subset(RA1_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA1_PBMC_seurat <- NormalizeData(RA1_PBMC_seurat, verbose = FALSE)
RA1_PBMC_seurat <- FindVariableFeatures(RA1_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA1_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA2_PBMC_seurat object
RA2_PBMC_seurat <- CreateSeuratObject(counts = RA2_PBMC.data, project = "RA2_PBMC", min.cells = 5, min.features = 500)
RA2_PBMC_seurat$samples <- "RA2_PBMC"
RA2_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA2_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA2_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA2_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA2_PBMC_seurat <- subset(RA2_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA2_PBMC_seurat <- NormalizeData(RA2_PBMC_seurat, verbose = FALSE)
RA2_PBMC_seurat <- FindVariableFeatures(RA2_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA2_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA3_PBMC_seurat object
RA3_PBMC_seurat <- CreateSeuratObject(counts = RA3_PBMC.data, project = "RA3_PBMC", min.cells = 5, min.features = 500)
RA3_PBMC_seurat$samples <- "RA3_PBMC"
RA3_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA3_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA3_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA3_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA3_PBMC_seurat <- subset(RA3_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA3_PBMC_seurat <- NormalizeData(RA3_PBMC_seurat, verbose = FALSE)
RA3_PBMC_seurat <- FindVariableFeatures(RA3_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA3_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA4_PBMC_seurat object
RA4_PBMC_seurat <- CreateSeuratObject(counts = RA4_PBMC.data, project = "RA4_PBMC", min.cells = 5, min.features = 500)
RA4_PBMC_seurat$samples <- "RA4_PBMC"
RA4_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA4_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA4_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA4_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA4_PBMC_seurat <- subset(RA4_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA4_PBMC_seurat <- NormalizeData(RA4_PBMC_seurat, verbose = FALSE)
RA4_PBMC_seurat <- FindVariableFeatures(RA4_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA4_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA5_PBMC_seurat object
RA5_PBMC_seurat <- CreateSeuratObject(counts = RA5_PBMC.data, project = "RA5_PBMC", min.cells = 5, min.features = 500)
RA5_PBMC_seurat$samples <- "RA5_PBMC"
RA5_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA5_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA5_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA5_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA5_PBMC_seurat <- subset(RA5_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA5_PBMC_seurat <- NormalizeData(RA5_PBMC_seurat, verbose = FALSE)
RA5_PBMC_seurat <- FindVariableFeatures(RA5_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA5_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA6_PBMC_seurat object
RA6_PBMC_seurat <- CreateSeuratObject(counts = RA6_PBMC.data, project = "RA6_PBMC", min.cells = 5, min.features = 500)
RA6_PBMC_seurat$samples <- "RA6_PBMC"
RA6_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA6_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA6_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA6_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA6_PBMC_seurat <- subset(RA6_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA6_PBMC_seurat <- NormalizeData(RA6_PBMC_seurat, verbose = FALSE)
RA6_PBMC_seurat <- FindVariableFeatures(RA6_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA6_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA8_PBMC_seurat object
RA8_PBMC_seurat <- CreateSeuratObject(counts = RA8_PBMC.data, project = "RA8_PBMC", min.cells = 5, min.features = 500)
RA8_PBMC_seurat$samples <- "RA8_PBMC"
RA8_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA8_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA8_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA8_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA8_PBMC_seurat <- subset(RA8_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA8_PBMC_seurat <- NormalizeData(RA8_PBMC_seurat, verbose = FALSE)
RA8_PBMC_seurat <- FindVariableFeatures(RA8_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA8_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA9_PBMC_seurat object
RA9_PBMC_seurat <- CreateSeuratObject(counts = RA9_PBMC.data, project = "RA9_PBMC", min.cells = 5, min.features = 500)
RA9_PBMC_seurat$samples <- "RA9_PBMC"
RA9_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA9_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA9_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA9_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA9_PBMC_seurat <- subset(RA9_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA9_PBMC_seurat <- NormalizeData(RA9_PBMC_seurat, verbose = FALSE)
RA9_PBMC_seurat <- FindVariableFeatures(RA9_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA9_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA10_PBMC_seurat object
RA10_PBMC_seurat <- CreateSeuratObject(counts = RA10_PBMC.data, project = "RA10_PBMC", min.cells = 5, min.features = 500)
RA10_PBMC_seurat$samples <- "RA10_PBMC"
RA10_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA10_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA10_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA10_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA10_PBMC_seurat <- subset(RA10_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA10_PBMC_seurat <- NormalizeData(RA10_PBMC_seurat, verbose = FALSE)
RA10_PBMC_seurat <- FindVariableFeatures(RA10_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA10_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA11_PBMC_seurat object
RA11_PBMC_seurat <- CreateSeuratObject(counts = RA11_PBMC.data, project = "RA11_PBMC", min.cells = 5, min.features = 500)
RA11_PBMC_seurat$samples <- "RA11_PBMC"
RA11_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA11_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA11_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA11_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA11_PBMC_seurat <- subset(RA11_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA11_PBMC_seurat <- NormalizeData(RA11_PBMC_seurat, verbose = FALSE)
RA11_PBMC_seurat <- FindVariableFeatures(RA11_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA11_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA12_PBMC_seurat object
RA12_PBMC_seurat <- CreateSeuratObject(counts = RA12_PBMC.data, project = "RA12_PBMC", min.cells = 5, min.features = 500)
RA12_PBMC_seurat$samples <- "RA12_PBMC"
RA12_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA12_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA12_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA12_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA12_PBMC_seurat <- subset(RA12_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA12_PBMC_seurat <- NormalizeData(RA12_PBMC_seurat, verbose = FALSE)
RA12_PBMC_seurat <- FindVariableFeatures(RA12_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA12_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA13_PBMC_seurat object
RA13_PBMC_seurat <- CreateSeuratObject(counts = RA13_PBMC.data, project = "RA13_PBMC", min.cells = 5, min.features = 500)
RA13_PBMC_seurat$samples <- "RA13_PBMC"
RA13_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA13_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA13_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA13_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA13_PBMC_seurat <- subset(RA13_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA13_PBMC_seurat <- NormalizeData(RA13_PBMC_seurat, verbose = FALSE)
RA13_PBMC_seurat <- FindVariableFeatures(RA13_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA13_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA14_PBMC_seurat object
RA14_PBMC_seurat <- CreateSeuratObject(counts = RA14_PBMC.data, project = "RA14_PBMC", min.cells = 5, min.features = 500)
RA14_PBMC_seurat$samples <- "RA14_PBMC"
RA14_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA14_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA14_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA14_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA14_PBMC_seurat <- subset(RA14_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA14_PBMC_seurat <- NormalizeData(RA14_PBMC_seurat, verbose = FALSE)
RA14_PBMC_seurat <- FindVariableFeatures(RA14_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA14_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA15_PBMC_seurat object
RA15_PBMC_seurat <- CreateSeuratObject(counts = RA15_PBMC.data, project = "RA15_PBMC", min.cells = 5, min.features = 500)
RA15_PBMC_seurat$samples <- "RA15_PBMC"
RA15_PBMC_seurat[["percent.mt"]] <- PercentageFeatureSet(RA15_PBMC_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA15_PBMC_seurat.pdf", width = 8, height = 6)
VlnPlot(RA15_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA15_PBMC_seurat <- subset(RA15_PBMC_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA15_PBMC_seurat <- NormalizeData(RA15_PBMC_seurat, verbose = FALSE)
RA15_PBMC_seurat <- FindVariableFeatures(RA15_PBMC_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA15_PBMC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#######################################################################################################
RA_PBMC_anchors <- FindIntegrationAnchors(object.list = list(HC1_PBMC_seurat, HC2_PBMC_seurat,HC3_PBMC_seurat,HC4_PBMC_seurat,RA1_PBMC_seurat, 
                                                             RA2_PBMC_seurat, RA3_PBMC_seurat, RA4_PBMC_seurat, RA5_PBMC_seurat, RA6_PBMC_seurat,
                                                             RA8_PBMC_seurat, RA9_PBMC_seurat, RA10_PBMC_seurat, RA11_PBMC_seurat, RA12_PBMC_seurat,
                                                             RA13_PBMC_seurat, RA14_PBMC_seurat, RA15_PBMC_seurat), 
                                          dims = 1:50,
                                          anchor.features = 3500)

save(RA_PBMC_anchors, file = "RA_PBMC_anchors.Rda")
RA_PBMC <- IntegrateData(anchorset = RA_PBMC_anchors, dims = 1:50)

save(RA_PBMC, file = "RA_PBMC_integrated_data.Rda")
RA_PBMC
DefaultAssay(RA_PBMC) <- "integrated"
#distinctCols <- distinctColorPalette(50)

########### Adding some features to the metadata---------------------------

current.cluster.ids <- c("HC1_PBMC", "HC2_PBMC", "HC3_PBMC","HC4_PBMC", "RA1_PBMC", "RA2_PBMC",
                         "RA3_PBMC", "RA4_PBMC", "RA5_PBMC", "RA6_PBMC", "RA8_PBMC", "RA9_PBMC", "RA10_PBMC", 
                         "RA11_PBMC", "RA12_PBMC", "RA13_PBMC", "RA14_PBMC", "RA15_PBMC")
new.cluster.ids <- c("HC","HC","HC","HC","Long","Long", "Long", "Long", "Long",
                     "Short", "Short", "Long", "Short", "Short", "Short", "Long", "Long", "Long")
RA_PBMC@meta.data$Duration <- plyr::mapvalues(x = RA_PBMC@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_PBMC@meta.data)
######
######
current.cluster.ids <- c("HC1_PBMC", "HC2_PBMC", "HC3_PBMC","HC4_PBMC", "RA1_PBMC", "RA2_PBMC",
                         "RA3_PBMC", "RA4_PBMC", "RA5_PBMC", "RA6_PBMC", "RA8_PBMC", "RA9_PBMC", "RA10_PBMC", 
                         "RA11_PBMC", "RA12_PBMC", "RA13_PBMC", "RA14_PBMC", "RA15_PBMC")
new.cluster.ids <- c("HC","HC","HC","HC","Positive","Negative", "Positive", "Positive", "Positive", 
                     "Positive","Positive", "Negative", "Negative", "Positive", "Negative", "Negative", "Negative", "Negative")
RA_PBMC@meta.data$ACPA <- plyr::mapvalues(x = RA_PBMC@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_PBMC@meta.data)
######
######
current.cluster.ids <- c("HC1_PBMC", "HC2_PBMC", "HC3_PBMC","HC4_PBMC", "RA1_PBMC", "RA2_PBMC",
                         "RA3_PBMC", "RA4_PBMC", "RA5_PBMC", "RA6_PBMC", "RA8_PBMC", "RA9_PBMC", "RA10_PBMC", 
                         "RA11_PBMC", "RA12_PBMC", "RA13_PBMC", "RA14_PBMC", "RA15_PBMC")
new.cluster.ids <- c("HC","HC","HC","HC","Pos_long","Neg_long", "Pos_long", "Pos_long", "Pos_long", "Pos_short",
                     "Pos_short", "Neg_long", "Neg_short", "Pos_short", "Neg_short", "Neg_long", "Neg_long", "Neg_long")
RA_PBMC@meta.data$Stage <- plyr::mapvalues(x = RA_PBMC@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_PBMC@meta.data)

#### Calculate the cell cycle score
cc.gene <- read.csv("./result/Seurat/cell_cycle_genes.csv",header=F)
s.genes <- cc.gene[1:43,1]
s.genes <- as.character(s.genes)
g2m.genes <- cc.gene[44:98,1]
g2m.genes <- as.character(g2m.genes)
# all(s.genes %in% rownames(RA_PBMC@data))
# all(g2m.genes %in% rownames(RA_PBMC@data))
# sum(s.genes %in% rownames(RA_PBMC@data))
# sum(g2m.genes %in% rownames(RA_PBMC@data))
s.genes <- intersect(s.genes, rownames(RA_PBMC))
g2m.genes <- intersect(g2m.genes, rownames(RA_PBMC))
RA_PBMC <- CellCycleScoring(object = RA_PBMC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(RA_PBMC@meta.data)
tail(RA_PBMC@meta.data)


# Run the standard workflow for visualization and clustering
RA_PBMC <- ScaleData(object = RA_PBMC, 
                     vars.to.regress = c("S.Score", "G2M.Score","percent.mt","nRNA_counts","nFeature_counts"), 
                     features = rownames(RA_PBMC))
RA_PBMC <- RunPCA(object = RA_PBMC,  
                  features = VariableFeatures(object = RA_PBMC), 
                  npcs = 50, 
                  verbose = FALSE)
RA_PBMC <- JackStraw(RA_PBMC, num.replicate = 100)
RA_PBMC <- ScoreJackStraw(RA_PBMC, dims = 1:20)
JackStrawPlot(RA_PBMC, dims = 1:20)
ElbowPlot(RA_PBMC)

RA_PBMC <- FindNeighbors(RA_PBMC, reduction = "pca", dims = 1:50)
RA_PBMC <- FindClusters(RA_PBMC, resolution = 0.8)

RA_PBMC <- RunUMAP(RA_PBMC, reduction = "pca", dims = 1:50)
RA_PBMC <- RunTSNE(RA_PBMC, reduction = "pca", dims = 1:50)

DimPlot(RA_PBMC,
        label = T, 
        label.size = 6, 
        cols = cell_type_cols)

save(RA_PBMC, file = "RA_PBMC_res0.8.Rda")

##############################################
##############################################
DefaultAssay(RA_PBMC) <- "RNA"
RA_PBMC.markers <- FindAllMarkers(RA_PBMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(RA_PBMC.markers, file = "./RA_PBMC_markers.Rdata")
x <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(x,file="RA_PBMC_top50_markers.csv")

x <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(x,file="RA_PBMC_top100_markers.csv")

x <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
write.csv(x,file="RA_PBMC_top200_markers.csv")

x <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
write.csv(x,file="RA_PBMC_top500_markers.csv")
#########################################################################################################
DefaultAssay(RA_PBMC) <- "RNA"
x <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
specific_genes <- unique(x$gene)

pdf(file = "Dotplot_RA_PBMC_features_genes.pdf", width = 16, height = 8)
DotPlot(object = RA_PBMC, 
        features = rev(specific_genes), 
        col.max = 2, 
        col.min = -1.5,
        cols = c("#436EEE","#FF6A6A","#FFA500","#FF34B3","#ADFF2F","#00F5FF")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#RColorBrewer::brewer.pal.info
#display.brewer.pal(9,"BrBG")
Dotplot_cols <- RColorBrewer::brewer.pal(8,"Set2")

pdf(file = "Dotplot_RA_PBMC_features_genes_split_by_ACPA.pdf", width = 18, height = 15)
DotPlot(object = RA_PBMC, 
        features = rev(specific_genes), 
        split.by = "ACPA",
        cols = Dotplot_cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


pdf(file = "Dotplot_RA_PBMC_features_genes_split_by_Stage.pdf", width = 18, height = 15)
DotPlot(object = RA_PBMC, 
        features = rev(specific_genes), 
        split.by = "Stage",
        cols = Dotplot_cols) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
#########################################################################################################
#########################################################################################################
###plot heatmap
DefaultAssay(RA_PBMC) <- "RNA"

top20 <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf(file="./heatmap_RA_PBMC_top20.pdf", width = 12, height = 10)
DoHeatmap(object = RA_PBMC,features = top20$gene, size = 3) +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                       mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()
DefaultAssay(RA_PBMC) <- "RNA"
################################################################################################
################################################################################################
pdf(file="./Umap_number.pdf", width = 8, height = 6)
DimPlot(object = RA_PBMC,
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        cols= cell_type_cols)

DimPlot(object = RA_PBMC,
        label = TRUE, 
        label.size = 6,
        pt.size = 0.1,
        split.by = "seurat_clusters",
        ncol = 6,
        cols= cell_type_cols)
dev.off()

pdf(file="./Umap_by_sample.pdf", width = 8, height = 6)
DimPlot(RA_PBMC, 
        label = F, 
        label.size = 4,
        pt.size = 0.5,
        group.by = "samples", 
        cols= cell_type_cols)
dev.off()

pdf(file="./Umap_by_sample_split.pdf", width = 12, height = 12)
DimPlot(object = RA_PBMC, label = F, 
        label.size = 4,
        pt.size = 0.1,
        group.by="samples",
        split.by = "samples",
        cols= cell_type_cols,
        ncol = 3)
dev.off()

pdf(file="./Umap_by_Phase.pdf", width = 8, height = 6)
DimPlot(object = RA_PBMC, reduction="umap",
        cols=c(brewer.pal(8, "Set2")),
        label = F, 
        label.size = 4,
        pt.size = 0.5, 
        group.by = "Phase", 
        order = rev(c("G1","S","G2M")))
dev.off()

pdf(file="./Umap_split_by_Phase.pdf", width = 10, height = 4)
DimPlot(object = RA_PBMC, reduction="umap",
        cols=c(brewer.pal(9, "Set1")),
        label = F, 
        label.size = 4,
        pt.size = 0.01, 
        group.by = "Phase",
        split.by = "Phase",
        ncol = 3,
        order = rev(c("G1","S","G2M")))
dev.off()

pdf(file="./Umap_by_ACPA_split_by_ACPA.pdf", width = 8, height = 3)
DimPlot(object = RA_PBMC, reduction="umap",
        cols=c(brewer.pal(8, "Set2")),
        label = F, 
        label.size = 4,
        pt.size = 0.1, 
        split.by = "ACPA",
        group.by = "ACPA")
dev.off()

pdf(file="./Umap_by_Stage.pdf", width = 8, height = 6)
DimPlot(object = RA_PBMC, reduction="umap",
        cols=c(brewer.pal(6, "Set2")),
        label = F, 
        label.size = 4,
        pt.size = 0.1, 
        group.by = "Stage")
dev.off()

pdf(file="./Umap_by_Stage_split_by_satge.pdf", width = 10, height = 2.5)
DimPlot(object = RA_PBMC, reduction="umap",
        cols=c(brewer.pal(8, "Set2")),
        label = F, 
        label.size = 4,
        pt.size = 0.1, 
        ncol = 5,
        split.by = "Stage",
        group.by = "Stage")
dev.off()

#####################################################################################
#####################################################################################
pdf(file = "Feature_specific_markers-1.pdf", width = 12, height = 8)
FeaturePlot(RA_PBMC,
            cols = c("grey80","red"),
            max.cutoff = 3, 
            ncol = 3,
            features = c("CD3E","CD4","CD8A","NKG7","GZMB","HLA-DRB5","CD68","HLA-DRB1","PPBP"))
dev.off()

pdf(file = "Feature_specific_markers-2.pdf", width = 12, height = 8)
FeaturePlot(RA_PBMC,
            cols = c("grey80","red"),
            max.cutoff = 3, 
            nncol = 3,
            features = c("CCR7","IRF8", "KLRC1", "GNLY","HBA2", "MX1", "ISG15", "SOX4","MS4A7"))
dev.off()


pdf(file="./Vlnplot_specific_markers-1.pdf", width = 12, height = 10)
VlnPlot(object = RA_PBMC, 
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0,
        features = c("CD3E","CD4","CD8A","NKG7","GZMB","HLA-DRB5","CD68","HLA-DRB1","PPBP"))
dev.off()

pdf(file="./Vlnplot_specific_markers-2.pdf", width = 12, height = 10)
VlnPlot(object = RA_PBMC, 
        cols = cell_type_cols,
        ncol = 3,
        pt.size = 0,
        features = c("CCR7","IRF8", "KLRC1", "GNLY","HBA2", "MX1", "ISG15", "SOX4","MS4A7"))
dev.off()
############################################################################
############################################################################
# current.cluster.ids <- c(0:20)
# new.cluster.ids <- c()
# RA_PBMC@meta.data$cell.type <- plyr::mapvalues(x = RA_PBMC@meta.data$seurat_clusters,
#                                                     from = current.cluster.ids,
#                                                     to = new.cluster.ids)
# head(RA_PBMC@meta.data)
# pdf(file="./Umap_by_cell.type.pdf", width = 10, height = 6)
# DimPlot(object = RA_PBMC,
#         label = T,
#         label.size = 4,
#         pt.size = 0.1,
#         group.by = "cell.type",
#         cols=cell_type_cols)
# 
# dev.off()

###########################################################################################
###########################################################################################
pdf("./Heatmap_ave_RA_PBMC_top80_cell_type.pdf", width = 8, height = 15)
#Idents(RA_PBMC) <- "seurat_clusters"
cluster.averages <- AverageExpression(RA_PBMC, return.seurat = T)
#my_cluster.averages <- cluster.averages[my_gene_set,]
DoHeatmap(cluster.averages, 
          features = unlist(TopFeatures(RA_PBMC[["pca"]], balanced = TRUE, nfeatures = 80, dim = 1)), 
          size = 4, draw.lines = FALSE, hjust = 0, group.bar = T) +
  theme(text = element_text(size = 15), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 10))+
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                       mid = "white",
                       high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill")

dev.off()
###########################################################################################
###########################################################################################
# Idents(object = RA_PBMC) <- "seurat_clusters"
av.exp <- AverageExpression(RA_PBMC)$RNA
Top_gene <- RA_PBMC.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
my_genes <- unique(Top_gene$gene)

av.exp_1 <- av.exp[my_genes,]
my_genes
cor.exp <- as.data.frame(cor(av.exp_1))

pdf(file="pheatmap_RA_PBMC-corr_top50_hcluster.pdf", width = 8, height = 8)
pheatmap::pheatmap(cor.exp,fontsize_row = 15,
                   fontsize_col = 15,
                   cluster_cols = T, 
                   angle_col = c(90),
                   cluster_rows = T,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=6,name="GnBu"))(60)) 
dev.off()

pdf(file="pheatmap_RA_PBMC-corr_top50.pdf")
pheatmap::pheatmap(cor.exp,fontsize_row = 15,
                   fontsize_col = 15,
                   cluster_cols = F, 
                   angle_col = c(90),
                   cluster_rows = F,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=6,name="YlGnBu"))(60)) 
dev.off()

############################################################################################
############################################################################################
## To check the percentage of cells in each cluster distributed into each organ.
Cluster_tissue_stat <- table( RA_PBMC@meta.data$samples, RA_PBMC@meta.data$seurat_clusters)
Cluster_tissue_stat
Cluster_tissue_stat_ordered <- apply(Cluster_tissue_stat, 1, function(x) x / sum(x))
Cluster_tissue_stat_ordered
#Cluster_tissue_stat_ordered <- t(Cluster_tissue_stat)

options(repr.plot.width=10, repr.plot.height=7)

# order_mat <- t(apply(Cluster_tissue_stat_ordered, 1, order))
# max_ind_vec <- c()
# 
# for(i in 1:nrow(order_mat)) {
#   tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
#   max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
#   }
# 
# max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
# 
# max_ind_vec <- c(max_ind_vec, setdiff(1:nrow(Cluster_tissue_stat),
#                                       max_ind_vec))
# Cluster_tissue_stat_ordered <- Cluster_tissue_stat_ordered[, row.names(Cluster_tissue_stat)[max_ind_vec]]

pdf(file="pheatmap_Fraction_of_cluster_per_sample.pdf", width = 6, height = 8)
pheatmap::pheatmap(Cluster_tissue_stat_ordered,
                   fontsize_row = 13,
                   fontsize_col = 10, 
                   cluster_cols = F, 
                   angle_col = c(315),
                   cluster_rows = F,
                   legend_labels = "Fraction of cell type",
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=5,name='BuGn'))(100)) 
dev.off()

pdf(file="pheatmap_Fraction_of_cluster_per_sample_hcluster.pdf", width = 6, height = 8)
pheatmap::pheatmap(Cluster_tissue_stat_ordered,
                   fontsize_row = 13,
                   fontsize_col = 10, 
                   cluster_cols = T, 
                   angle_col = c(315),
                   cluster_rows = T,
                   legend_labels = "Fraction of cell type",
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=5,name='BuGn'))(100)) 
dev.off()

############################################################################################
############################################################################################

## To check the percentage of cells in each cluster distributed into each organ.
Cluster_tissue_stat <- table(RA_PBMC@meta.data$seurat_clusters, RA_PBMC@meta.data$ACPA)
Cluster_tissue_stat
Cluster_tissue_stat_ordered <- apply(Cluster_tissue_stat, 1, function(x) x / sum(x))
Cluster_tissue_stat_ordered
#Cluster_tissue_stat_ordered <- t(Cluster_tissue_stat)

options(repr.plot.width=10, repr.plot.height=7)

# order_mat <- t(apply(Cluster_tissue_stat_ordered, 1, order))
# max_ind_vec <- c()
# 
# for(i in 1:nrow(order_mat)) {
#   tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
#   max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
#   }
# 
# max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
# 
# max_ind_vec <- c(max_ind_vec, setdiff(1:nrow(Cluster_tissue_stat),
#                                       max_ind_vec))
# Cluster_tissue_stat_ordered <- Cluster_tissue_stat_ordered[, row.names(Cluster_tissue_stat)[max_ind_vec]]

pdf(file="pheatmap_Fraction_of_cluster_by_ACPA.pdf", width = 8, height = 3)
pheatmap::pheatmap(Cluster_tissue_stat_ordered,
                   fontsize_row = 13,
                   fontsize_col = 10, 
                   cluster_cols = F, 
                   angle_col = c(0),
                   cluster_rows = F,
                   legend_labels = "Fraction of cell type",
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=5,name='BuGn'))(100)) 
dev.off()


######################
##SM
#/User/bin/env Rscript
## Seurat standard workplow
library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)
library(umap)
library(RColorBrewer)
library(cowplot)
library(randomcoloR)
library(corrplot)
library(DoubletFinder)

cell_type_cols <- c(brewer.pal(9, "Set1"),
                    brewer.pal(8,"Dark2"),
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4",
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030",
                    "#7CFC00","#000000","#708090")

RA2_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA2-SM/GRCh38")
RA3_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA3-SM/GRCh38")
RA4_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA4-SM/GRCh38/")
RA5_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA5-SM/GRCh38/")
RA6_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA6-SM/GRCh38/")
RA8_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA08-SM/GRCh38")
RA9_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA09-SM/GRCh38")
RA10_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA10-SM/GRCh38")
RA11_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA11-SM/GRCh38")
RA12_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA12-SM/GRCh38")
RA13_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA13-SM/GRCh38/")
RA14_SM.data <- Read10X(data.dir = "D:/RA/samples/SM/RA14-SM/GRCh38/")

# Singe Seurat object built -----------------------------------------------
# Set up RA2_SM_seurat object
RA2_SM_seurat <- CreateSeuratObject(counts = RA2_SM.data, project = "RA2_SM", min.cells = 5, min.features = 500)
RA2_SM_seurat$samples <- "RA2_SM"
RA2_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA2_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA2_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA2_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA2_SM_seurat <- subset(RA2_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA2_SM_seurat <- NormalizeData(RA2_SM_seurat, verbose = FALSE)
RA2_SM_seurat <- FindVariableFeatures(RA2_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA2_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA3_SM_seurat object
RA3_SM_seurat <- CreateSeuratObject(counts = RA3_SM.data, project = "RA3_SM", min.cells = 5, min.features = 500)
RA3_SM_seurat$samples <- "RA3_SM"
RA3_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA3_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA3_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA3_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA3_SM_seurat <- subset(RA3_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA3_SM_seurat <- NormalizeData(RA3_SM_seurat, verbose = FALSE)
RA3_SM_seurat <- FindVariableFeatures(RA3_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA3_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA4_SM_seurat object
RA4_SM_seurat <- CreateSeuratObject(counts = RA4_SM.data, project = "RA4_SM", min.cells = 5, min.features = 500)
RA4_SM_seurat$samples <- "RA4_SM"
RA4_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA4_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA4_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA4_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA4_SM_seurat <- subset(RA4_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA4_SM_seurat <- NormalizeData(RA4_SM_seurat, verbose = FALSE)
RA4_SM_seurat <- FindVariableFeatures(RA4_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA4_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA5_SM_seurat object
RA5_SM_seurat <- CreateSeuratObject(counts = RA5_SM.data, project = "RA5_SM", min.cells = 5, min.features = 500)
RA5_SM_seurat$samples <- "RA5_SM"
RA5_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA5_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA5_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA5_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA5_SM_seurat <- subset(RA5_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA5_SM_seurat <- NormalizeData(RA5_SM_seurat, verbose = FALSE)
RA5_SM_seurat <- FindVariableFeatures(RA5_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA5_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA6_SM_seurat object
RA6_SM_seurat <- CreateSeuratObject(counts = RA6_SM.data, project = "RA6_SM", min.cells = 5, min.features = 500)
RA6_SM_seurat$samples <- "RA6_SM"
RA6_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA6_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA6_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA6_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA6_SM_seurat <- subset(RA6_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA6_SM_seurat <- NormalizeData(RA6_SM_seurat, verbose = FALSE)
RA6_SM_seurat <- FindVariableFeatures(RA6_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA6_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA8_SM_seurat object
RA8_SM_seurat <- CreateSeuratObject(counts = RA8_SM.data, project = "RA8_SM", min.cells = 5, min.features = 500)
RA8_SM_seurat$samples <- "RA8_SM"
RA8_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA8_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA8_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA8_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA8_SM_seurat <- subset(RA8_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA8_SM_seurat <- NormalizeData(RA8_SM_seurat, verbose = FALSE)
RA8_SM_seurat <- FindVariableFeatures(RA8_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA8_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA9_SM_seurat object
RA9_SM_seurat <- CreateSeuratObject(counts = RA9_SM.data, project = "RA9_SM", min.cells = 5, min.features = 500)
RA9_SM_seurat$samples <- "RA9_SM"
RA9_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA9_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA9_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA9_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA9_SM_seurat <- subset(RA9_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA9_SM_seurat <- NormalizeData(RA9_SM_seurat, verbose = FALSE)
RA9_SM_seurat <- FindVariableFeatures(RA9_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA9_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA10_SM_seurat object
RA10_SM_seurat <- CreateSeuratObject(counts = RA10_SM.data, project = "RA10_SM", min.cells = 5, min.features = 500)
RA10_SM_seurat$samples <- "RA10_SM"
RA10_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA10_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA10_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA10_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA10_SM_seurat <- subset(RA10_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA10_SM_seurat <- NormalizeData(RA10_SM_seurat, verbose = FALSE)
RA10_SM_seurat <- FindVariableFeatures(RA10_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA10_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA11_SM_seurat object
RA11_SM_seurat <- CreateSeuratObject(counts = RA11_SM.data, project = "RA11_SM", min.cells = 5, min.features = 500)
RA11_SM_seurat$samples <- "RA11_SM"
RA11_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA11_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA11_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA11_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA11_SM_seurat <- subset(RA11_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA11_SM_seurat <- NormalizeData(RA11_SM_seurat, verbose = FALSE)
RA11_SM_seurat <- FindVariableFeatures(RA11_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA11_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA12_SM_seurat object
RA12_SM_seurat <- CreateSeuratObject(counts = RA12_SM.data, project = "RA12_SM", min.cells = 5, min.features = 500)
RA12_SM_seurat$samples <- "RA12_SM"
RA12_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA12_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA12_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA12_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA12_SM_seurat <- subset(RA12_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA12_SM_seurat <- NormalizeData(RA12_SM_seurat, verbose = FALSE)
RA12_SM_seurat <- FindVariableFeatures(RA12_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA12_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA13_SM_seurat object
RA13_SM_seurat <- CreateSeuratObject(counts = RA13_SM.data, project = "RA13_SM", min.cells = 5, min.features = 500)
RA13_SM_seurat$samples <- "RA13_SM"
RA13_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA13_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA13_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA13_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA13_SM_seurat <- subset(RA13_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA13_SM_seurat <- NormalizeData(RA13_SM_seurat, verbose = FALSE)
RA13_SM_seurat <- FindVariableFeatures(RA13_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA13_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Set up RA14_SM_seurat object
RA14_SM_seurat <- CreateSeuratObject(counts = RA14_SM.data, project = "RA14_SM", min.cells = 5, min.features = 500)
RA14_SM_seurat$samples <- "RA14_SM"
RA14_SM_seurat[["percent.mt"]] <- PercentageFeatureSet(RA14_SM_seurat, pattern = "^MT-")
pdf(file="Vlnplot_RA14_SM_seurat.pdf", width = 8, height = 6)
VlnPlot(RA14_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RA14_SM_seurat <- subset(RA14_SM_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5 & nCount_RNA < 20000)
RA14_SM_seurat <- NormalizeData(RA14_SM_seurat, verbose = FALSE)
RA14_SM_seurat <- FindVariableFeatures(RA14_SM_seurat, selection.method = "vst", nfeatures = 3500)
VlnPlot(RA14_SM_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#######################################################################################################
RA_SM_anchors <- FindIntegrationAnchors(object.list = list( RA2_SM_seurat, RA3_SM_seurat,RA5_SM_seurat, RA6_SM_seurat,
                                                            RA8_SM_seurat, RA9_SM_seurat, RA10_SM_seurat, RA11_SM_seurat, RA12_SM_seurat,
                                                            RA14_SM_seurat), 
                                        dims = 1:50,
                                        anchor.features = 3500)

save(RA_SM_anchors, file = "RA_SM_anchors.Rda")
RA_SM <- IntegrateData(anchorset = RA_SM_anchors, dims = 1:50)
RA_SM
DefaultAssay(RA_SM) <- "integrated"
#distinctCols <- distinctColorPalette(50)

########### Adding some features to the metadata---------------------------

current.cluster.ids <- c("RA2_SM",
                         "RA3_SM", "RA4_SM", "RA5_SM", "RA6_SM", "RA8_SM", "RA9_SM", "RA10_SM", 
                         "RA11_SM", "RA12_SM", "RA13_SM", "RA14_SM")
new.cluster.ids <- c("Long", "Long", "Long", "Long",
                     "Short", "Short", "Long", "Short", "Short", "Short", "Long", "Long")
RA_SM@meta.data$Duration <- plyr::mapvalues(x = RA_SM@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_SM@meta.data)
######
######
current.cluster.ids <- c("RA2_SM",
                         "RA3_SM", "RA4_SM", "RA5_SM", "RA6_SM", "RA8_SM", "RA9_SM", "RA10_SM", 
                         "RA11_SM", "RA12_SM", "RA13_SM", "RA14_SM")
new.cluster.ids <- c("Positive", "Positive", "Positive", "Positive", 
                     "Positive","Positive", "Negative", "Negative", "Positive", "Negative", "Negative", "Negative")
RA_SM@meta.data$ACPA <- plyr::mapvalues(x = RA_SM@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_SM@meta.data)
######
######
current.cluster.ids <- c("RA2_SM","RA3_SM", "RA4_SM", "RA5_SM", "RA6_SM", "RA8_SM", "RA9_SM", "RA10_SM", 
                         "RA11_SM", "RA12_SM", "RA13_SM", "RA14_SM")
new.cluster.ids <- c("Pos_long", "Pos_long", "Pos_long", "Pos_long", "Pos_short",
                     "Pos_short", "Neg_long", "Neg_short", "Pos_short", "Neg_short", "Neg_long", "Neg_long")
RA_SM@meta.data$Stage <- plyr::mapvalues(x = RA_SM@meta.data$samples, from = current.cluster.ids, to = new.cluster.ids)
head(RA_SM@meta.data)

#### Calculate the cell cycle score
cc.gene <- read.csv("C:/Users/aking/Desktop/MF/cell_cycle_genes.csv",header=F)
# cc.gene <- read.csv("/Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/Seurat/cell_cycle_genes.csv",header=F)
s.genes <- cc.gene[1:43,1]
s.genes <- as.character(s.genes)
g2m.genes <- cc.gene[44:98,1]
g2m.genes <- as.character(g2m.genes)
# all(s.genes %in% rownames(RA_SM@data))
# all(g2m.genes %in% rownames(RA_SM@data))
# sum(s.genes %in% rownames(RA_SM@data))
# sum(g2m.genes %in% rownames(RA_SM@data))
s.genes <- intersect(s.genes, rownames(RA_SM))
g2m.genes <- intersect(g2m.genes, rownames(RA_SM))
RA_SM <- CellCycleScoring(object = RA_SM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(RA_SM@meta.data)
tail(RA_SM@meta.data)


# Run the standard workflow for visualization and clustering
RA_SM <- ScaleData(object = RA_SM, 
                   vars.to.regress = c("S.Score", "G2M.Score","percent.mt","nRNA_counts","nFeature_counts"), 
                   features = rownames(RA_SM))
RA_SM <- RunPCA(object = RA_SM,  
                features = VariableFeatures(object = RA_SM), 
                npcs = 50, 
                verbose = FALSE)
RA_SM <- JackStraw(RA_SM, num.replicate = 100)
RA_SM <- ScoreJackStraw(RA_SM, dims = 1:20)
JackStrawPlot(RA_SM, dims = 1:20)
ElbowPlot(RA_SM)

RA_SM <- FindNeighbors(RA_SM, reduction = "pca", dims = 1:50)
# RA_SM <- FindClusters(RA_SM, resolution = 0.8)
# 
# RA_SM <- RunUMAP(RA_SM, reduction = "pca", dims = 1:50)
# RA_SM <- RunTSNE(RA_SM, reduction = "pca", dims = 1:50)
RA_SMnew <- FindClusters(RA_SM, resolution = 0.5)

RA_SMnew <- RunUMAP(RA_Snew, reduction = "pca", dims = 1:50)
RA_SMnew <- RunTSNE(RA_SMnew, reduction = "pca", dims = 1:50)

DimPlot(RA_SMnew,
        label = T, 
        label.size = 6, 
        cols = cell_type_cols)

#save(RA_SM, file = "RA_SM_res0.8.Rda")
save(RA_SMnew, file = "RA_SMnew-30-0.5.RData")
##############################################
##############################################
DefaultAssay(RA_SM) <- "RNA"
RA_SM.markers <- FindAllMarkers(RA_SM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(RA_SM.markers, file = "./RA_SM_markers.Rdata")
x <- RA_SM.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(x,file="RA_SM_top50_markers.csv")

x <- RA_SM.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(x,file="RA_SM_top100_markers.csv")

x <- RA_SM.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
write.csv(x,file="RA_SM_top200_markers.csv")

x <- RA_SM.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
write.csv(x,file="RA_SM_top500_markers.csv")
#########################################################################################################
