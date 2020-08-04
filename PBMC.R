library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
load("RA_PBMC-labeld-0225.RData")
celltype1 <- c("Naïve T","Exhausted T","IFN-act T","Memory T",
               "Naïve CD8T","GZMKhi T","CXCL13hi T", "Cytotoxic T",
               "NKT","NK",
               "Naïve B","Mem-unsw B","Mem-sw B","plasma B","IRF8low B",
               "CD14+ M","IL1Bhi CD14+M","IRF7+ IFN-act CD14+M",
               "CD16+ M","IL1Bhi CD16+M","IRF8hi DC","cDC","pDC",
               "CD34+ cells","Dividing Cells","Red cells","platelet")
celltype2 <- c("Naïve T","Exhausted T","IFN-act T","Memory T",
               "Naïve CD8T","GZMKhi T","CXCR5+ T", "Cytotoxic T",
               "NKT","NK",
               "Naïve B","Mem-unsw B","Mem-sw B","plasma B","CD20+CD14+ cells",
               "CD14+ M","HLA-DRB5+ CD14+M","IFN-act CD14+M",
               "CD16+ M","HLA-DRB5+ CD16+M","S100Bhi DC","cDC","pDC",
               "CD34+ cells","Dividing Cells","Red cells","platelet")
RA_PBMC@meta.data$new_celltype <- plyr::mapvalues(x=RA_PBMC@meta.data$new_celltype,
                                                  from=celltype1,to=celltype2)
RA_PBMC@meta.data$numbers <- plyr::mapvalues(x=RA_PBMC@meta.data$new_celltype,from = celltype2,
                                             to=c(1:27))
table(RA_PBMC@meta.data$new_celltype)
color27 <- c("#1e8bc3","#1f3a93","#75DDDD","#118ab2","#9faffd",
             "#5333ed","#c5eff7","#55e9ff","#446cb3","#59abe3",
             "#d44c00","#ef9746","#f5ba08","#eea6b9","#774936",
             "#849324","#5a622d","#c2a31a","#7ebc89","#1a936f","#52796f","#c6dabf","#127475",
             "#ffd670","#ffc43d","#fff3b0","#d8973c")

pdf("PBMC-0524/PBMC-umap.pdf",width = 8,height = 6)
DimPlot(RA_PBMC,reduction = "umap",group.by = "new_celltype",cols = color27)
DimPlot(RA_PBMC,reduction = "umap",group.by="numbers",cols = color27,label = T,label.size = 5)+
  ggmin::theme_powerpoint()+
  theme(panel.border = element_rect(fill = NA, colour = "#4e5569"))+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())+
  theme(legend.position = "none")
dev.off()
####pheatmap
Idents(RA_PBMC)<-"new_celltype"
pbmc_B <- subset(RA_PBMC,idents = c("Naïve B","Mem-unsw B","Mem-sw B","plasma B"))
levels(pbmc_B)
pbmc_B@meta.data$new_celltype <- factor(pbmc_B@meta.data$new_celltype,
                                        levels = c("Naïve B","Mem-unsw B","Mem-sw B","plasma B"))
Idents(pbmc_B) <- "new_celltype"
features <- c("IGHD","CD27","MS4A1","GNLY","MZB1","JCHAIN","IGHM","IGHG1","IGHG3","IGHA1")
DefaultAssay(pbmc_B) <- "RNA"
intersect(rownames(pbmc_B),features)
b.expr <- AverageExpression(object = pbmc_B,features = features,assays = "RNA")
pdf("PBMC-0524/pbmc_B_pheatmap.pdf",width = 5,height = 5)
pheatmap(b.expr$RNA,scale="row",cluster_rows = F, cluster_cols = F,angle_col = 315,fontsize=12)
dev.off()
##
pbmc_T <- subset(RA_PBMC,idents=c("Naïve T","Exhausted T","IFN-act T","Memory T",
                                  "Naïve CD8T","GZMKhi T","CXCR5+ T", "Cytotoxic T"))
levels(pbmc_T)
pbmc_T@meta.data$new_celltype <- factor(pbmc_T@meta.data$new_celltype,
                                        levels = c("Naïve T","Exhausted T","IFN-act T","Memory T",
                                                   "Naïve CD8T","GZMKhi T","CXCR5+ T", "Cytotoxic T"))
Idents(pbmc_T) <- "new_celltype"
features <- c("CD8A","NKG7","GNLY","GZMB","GZMK",
              "TBX21","IFNG","IFI6","ISG15","STAT1",
              "CCR7","SELL","IL7R","CD44","GPR183",
              "PDCD1","TIGIT","LAG3","ICOS","CXCR5")
DefaultAssay(pbmc_T) <- "RNA"
intersect(rownames(pbmc_T),features)
b.expr <- AverageExpression(object = pbmc_T,features = features,assays = "RNA")
pdf("PBMC-0524/pbmc_T_pheatmap.pdf",width = 5,height = 10)
pheatmap(b.expr$RNA,scale="row",cluster_rows = F, cluster_cols = F,angle_col = 315,fontsize=12)
dev.off()
##
pbmc_mono <- subset(RA_PBMC,idents=c("CD14+ M","HLA-DRB5+ CD14+M","IFN-act CD14+M",
                                     "CD16+ M","HLA-DRB5+ CD16+M","S100Bhi DC","cDC","pDC"))
levels(pbmc_mono)
pbmc_mono@meta.data$new_celltype <- factor(pbmc_mono@meta.data$new_celltype,
                                           levels = c("CD14+ M","HLA-DRB5+ CD14+M","IFN-act CD14+M",
                                                      "CD16+ M","HLA-DRB5+ CD16+M","S100Bhi DC","cDC","pDC"))
Idents(pbmc_mono) <- "new_celltype"
features <- c("CD14","FCGR3A","FCER1A","CD1C","IFI44","STAT1","IL1B",
              "ISG15","IFI6","HLA-DRB5","IRF8","IL3RA","S100B")
DefaultAssay(pbmc_mono) <- "RNA"
intersect(rownames(pbmc_mono),features)
b.expr <- AverageExpression(object = pbmc_mono,features = features,assays = "RNA")
pdf("PBMC-0524/pbmc_mono_pheatmap.pdf",width = 5,height = 8)
pheatmap(b.expr$RNA,scale="row",cluster_rows = F, cluster_cols = F,angle_col = 315,fontsize=12)
dev.off()
###vlnplot
features <- c("VPREB3","IGHD","CD79A","CD79B","MS4A1","JCHAIN","IGHA1","IGHG1","IGHM","CD14","CD1C","GPR183","CD22","CD19")
features<- intersect(features,rownames(RA_PBMC))
single_violin=function(seu.obj,feature,colors){
  VlnPlot(seu.obj,features = feature,group.by = "new_celltype",
          pt.size = 0,ncol = 1,cols = colors)+
    coord_flip()+
    theme(plot.title = element_text(size=8,hjust = 0.5))+
    theme(panel.background=element_rect(fill="white",color="white"))+
    scale_y_discrete(expand = c(0, 0))+
    theme(axis.text.y = element_blank())+
    theme(axis.title.y = element_blank(),axis.title.x = element_blank())+
    theme(axis.ticks.y = element_blank())+
    theme(legend.position = "NULL")+
    theme(plot.margin = unit(c(0,0,0,0),"mm"))+
    theme(panel.spacing = unit(c(0,0,0,0),"mm"))+
    theme(strip.placement = "inside")
}
###
PBMC_metadata <- RA_PBMC@meta.data
save(PBMC_metadata,file = "PBMC_metadata.RData")
#####
