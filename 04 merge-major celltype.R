## 
### ---------------
###
### Create: Yi Liu
### Date: 2019/12/02
### Update Log: 2020/08/04  
### This script aims to integrate major cell types from PBMC and SM for further sub-clustering 
###
### ---------------

###MonoDC
#load("./RData/pbmc_monodc.RData")
#load("./RData/sm_monodc.RData")
library(Seurat)
DefaultAssay(pbmc_monodc) <- "RNA"
DefaultAssay(sm_monodc) <- "RNA"
RA_Monodc_anchors <- FindIntegrationAnchors(object.list = list(pbmc_monodc,sm_monodc),
                                            dims = 1:50,
                                            anchor.features = 3500)
#save(RA_Monodc_anchors, file = "RA_Monodc_anchors-noSF.Rda")
#load("RA_Monodc_anchors-noSF.Rda")
RA_Monodc <- IntegrateData(anchorset = RA_Monodc_anchors, dims = 1:30)
RA_Monodc
DefaultAssay(RA_Monodc) <- "integrated"
save(RA_Monodc,file="RA_Monodc_just_integrated.Rda")
######################################
s.genes <- intersect(cc.genes$s.genes, rownames(RA_Monodc))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(RA_Monodc))
RA_Monodc <- CellCycleScoring(object = RA_Monodc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RA_Monodc <- ScaleData(object = RA_Monodc, 
                       vars.to.regress = c("S.Score", "G2M.Score","percent.mt","nRNA_counts","nFeature_counts"), 
                       features = rownames(RA_Monodc))
RA_Monodc <- RunPCA(object = RA_Monodc,  
                    features = VariableFeatures(object = RA_Monodc), 
                    npcs = 50, 
                    verbose = FALSE)
ElbowPlot(RA_Monodc)

RA_Monodc <- FindNeighbors(RA_Monodc, reduction = "pca", dims = 1:30)
RA_Monodc <- FindClusters(RA_Monodc, resolution = 0.8)

RA_Monodc <- RunUMAP(RA_Monodc, reduction = "pca", dims = 1:30)
RA_Monodc <- RunTSNE(RA_Monodc, reduction = "pca", dims = 1:30)
save(RA_Monodc, file = "RA_Monodc_res0.8_dim30-noSF.Rda")
pdf("pics3/monodc-dimplots-noSF.pdf",width = 10,height = 7)
DimPlot(RA_Monodc,label = T,label.size = 6)
DimPlot(RA_Monodc,label = T,reduction = "tsne",label.size = 6)
dev.off()
###T
load("./RData/pbmc_T.RData")
load("./RData/sm_t.RData")
library(Seurat)
DefaultAssay(sm_t) <- "RNA"
DefaultAssay(pbmc_T) <- "RNA"

RA_T_anchors <- FindIntegrationAnchors(object.list = list(pbmc_T,sm_t),
                                       dims = 1:50,
                                       anchor.features = 3500)
#save(RA_T_anchors, file = "RA_T_anchors-noSF.Rda")
#load("RA_T_anchors-noSF.Rda")
RA_T <- IntegrateData(anchorset = RA_T_anchors, dims = 1:50)
RA_T
DefaultAssay(RA_T) <- "integrated"
save(RA_T,file="RA_T_justmerged.RData")
print("...")
######################################
s.genes <- intersect(cc.genes$s.genes, rownames(RA_T))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(RA_T))
RA_T <- CellCycleScoring(object = RA_T, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RA_T <- ScaleData(object = RA_T, 
                  vars.to.regress = c("S.Score", "G2M.Score","percent.mt","nRNA_counts","nFeature_counts"), 
                  features = rownames(RA_T))
RA_T <- RunPCA(object = RA_T,  
               features = VariableFeatures(object = RA_T), 
               npcs = 50, 
               verbose = FALSE)
ElbowPlot(RA_T)

RA_T <- FindNeighbors(RA_T, reduction = "pca", dims = 1:30)
RA_T <- FindClusters(RA_T, resolution = 0.8)

RA_T <- RunUMAP(RA_T, reduction = "pca", dims = 1:30)
RA_T <- RunTSNE(RA_T, reduction = "pca", dims = 1:30)
save(RA_T, file = "RA_T_res0.8_dim30-noSF.Rda")
pdf("pics3/T-dimplots-noSF.pdf",width = 10,height = 7)
DimPlot(RA_T,label = T,label.size = 6)
DimPlot(RA_T,label = T,reduction="tsne",label.size = 6)
dev.off()
###B
load("./RData/pbmc_B.RData")
load("./RData/sm_b.RData")
DefaultAssay(pbmc_B) <- "RNA"
DefaultAssay(sm_b) <- "RNA"
RA_B_anchors <- FindIntegrationAnchors(object.list = list(pbmc_B,sm_b),
                                       dims = 1:50,
                                       anchor.features = 3500)
#save(RA_B_anchors, file = "RA_B_anchors-v2.Rda")
RA_B <- IntegrateData(anchorset = RA_B_anchors, dims = 1:30)
RA_B
DefaultAssay(RA_B) <- "integrated"
######################################
s.genes <- intersect(cc.genes$s.genes, rownames(RA_B))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(RA_B))
RA_B <- CellCycleScoring(object = RA_B, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RA_B <- ScaleData(object = RA_B, 
                  vars.to.regress = c("S.Score", "G2M.Score","percent.mt","nRNA_counts","nFeature_counts"), 
                  features = rownames(RA_B))
RA_B <- RunPCA(object = RA_B,  
               features = VariableFeatures(object = RA_B), 
               npcs = 50, 
               verbose = FALSE)
ElbowPlot(RA_B)

RA_B <- FindNeighbors(RA_B, reduction = "pca", dims = 1:30)
RA_B <- FindClusters(RA_B, resolution = 0.8)

RA_B <- RunUMAP(RA_B, reduction = "pca", dims = 1:30)
RA_B <- RunTSNE(RA_B, reduction = "pca", dims = 1:30)

save(RA_B, file = "RA_B_res0.8_dim30.Rda")
pdf("pics3/RA_B-umap.pdf",width = 10,height = 7)
DimPlot(RA_B,label = T,label.size = 6)
dev.off()