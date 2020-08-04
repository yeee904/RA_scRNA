## 
### ---------------
###
### Create: Yi Liu
### Date: 2019/11/30
### Update Log: 2020/08/04  
### This script aims to extract major cell types from PBMC and SM.
###
### ---------------

###setwd("/Bailab5/liuyi/RA-1210")
load("RA_PBMC-30-1.0.RData")
library(Seurat)
library(ggplot2)
RA_PBMC <- AddMetaData(RA_PBMC,metadata = Idents(RA_PBMC),col.name = "celltype")
RA_PBMC@meta.data$celltype <- plyr::mapvalues(x=RA_PBMC@meta.data$celltype,
                                              from = "IRF8neg pre-B",to="IRF8neg naïveB")
levels(RA_PBMC@meta.data$celltype)
RA_PBMC@meta.data$celltype <- factor(RA_PBMC@meta.data$celltype,
                                     levels = c("Naïve CD4T","Memory CD4T","IRF7+ memory CD4T","Effector CD4T",
                                                "Naïve CD8T","GZMK+ CD8T","IFN-act GZMK+ CD8T","CTL",
                                                "NKT","NK",
                                                "Naïve B","Memory-unswitched B","Memory-switched B","plasma B","IRF8neg naïveB",
                                                "CD14+ M","HLA-DRB5+ CD14+M","IRF7+ IFN-act CD14+M",
                                                "CD16+ M","HLA-DRB5+ CD16+M","IRFhi DC","cDC","pDC",
                                                "undeff-HSC","Dividing Cells","Red cells","platelet"))
RA_PBMC@meta.data$numbers <- plyr::mapvalues(x=RA_PBMC@meta.data$celltype,from = c("Naïve CD4T","Memory CD4T","IRF7+ memory CD4T","Effector CD4T",
                                                                                   "Naïve CD8T","GZMK+ CD8T","IFN-act GZMK+ CD8T","CTL",
                                                                                   "NKT","NK",
                                                                                   "Naïve B","Memory-unswitched B","Memory-switched B","plasma B","IRF8neg naïveB",
                                                                                   "CD14+ M","HLA-DRB5+ CD14+M","IRF7+ IFN-act CD14+M",
                                                                                   "CD16+ M","HLA-DRB5+ CD16+M","IRFhi DC","cDC","pDC",
                                                                                   "undeff-HSC","Dividing Cells","Red cells","platelet"),
                                             to=c(1:27))
save(RA_PBMC,file="RA_PBMC-labeled.RData")
Idents(RA_PBMC)<-"celltype"
pbmc_B <- subset(RA_PBMC,idents = c("Naïve B","Memory-unswitched B","Memory-switched B","plasma B","IRF8neg naïveB"))
save(pbmc_B,file="pbmc_B.RData")
pbmc_T <- subset(RA_PBMC,idents=c("Naïve CD4T","Memory CD4T","IRF7+ memory CD4T","Effector CD4T",
                                  "Naïve CD8T","GZMK+ CD8T","IFN-act GZMK+ CD8T","CTL",
                                  "NKT","NK"))
pbmc_monodc <- subset(RA_PBMC,idents=c("CD14+ M","HLA-DRB5+ CD14+M","IRF7+ IFN-act CD14+M",
                                       "CD16+ M","HLA-DRB5+ CD16+M","IRFhi DC","cDC","pDC"))
save(pbmc_T,file="pbmc_T.RData")
save(pbmc_monodc,file = "pbmc_monodc.RData")
##
load("SMnew-V2-30-1.0.RData")
levels(Idents(RA_SMnew))
celltypes <- c("CCL+ Mk","CCL18+MMP3+ Mk","CCL18+MMP3+ Mk","HLA-DRB5+ Mk","TEM",
               "IGHGhi plasmaB","CCL+ Mk","CCL+ Mk","CCL+ Mk","TEM",
               "Dividing cells","cDC","IL1Bhi Mk","Exhausted T","Memory B",
               "pDC","S100Bhi DC","LAMP3+ DC","Cytotoxic T","IFN-act Mk",
               "CD34+ cells","mast cells","plasma B")
RA_SMnew@meta.data$cell.type <- plyr::mapvalues(RA_SMnew@meta.data$seurat_clusters,
                                                from = c(0:22),to=celltypes)
celltype <- c("TEM","Cytotoxic T","Exhausted T",
              "Memory B","plasma B","IGHGhi plasmaB",
              "cDC","S100Bhi DC","LAMP3+ DC","pDC",
              "CCL+ Mk","CCL18+MMP3+ Mk","HLA-DRB5+ Mk","IFN-act Mk","IL1Bhi Mk",
              "CD34+ cells","mast cells","Dividing cells")
RA_SMnew@meta.data$cell.type <- factor(RA_SMnew@meta.data$cell.type,levels = celltype)
head(RA_SMnew@meta.data)
save(RA_SMnew,file="RA_SMnew-labeled.RData")
Idents(RA_SMnew) <- "cell.type"
bcell <- c("Memory B","plasma B","IGHGhi plasmaB")
tcell <- c("TEM","Cytotoxic T","Exhausted T")
dc <- c("cDC","S100Bhi DC","LAMP3+ DC","pDC")
mk <- c("CCL+ Mk","CCL18+MMP3+ Mk","HLA-DRB5+ Mk","IFN-act Mk","IL1Bhi Mk")
sm_t <- subset(RA_SMnew,idents=tcell)
sm_b <- subset(RA_SMnew,idents=bcell)
sm_monodc <- subset(RA_SMnew,idents=c(mk,dc))
save(sm_t,file="sm_t.RData")
save(sm_b,file="sm_b.RData")
save(sm_monodc,file="sm_monodc.RData")