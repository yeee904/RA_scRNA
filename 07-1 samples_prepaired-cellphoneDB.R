library(Seurat)
library(dplyr)

### Read in expression data of interacting cells
load("D:/我的课题/RA/0321上传/B_labeled_0321.RData")
load("D:/我的课题/RA/0321上传/Mk-labeled-0207.RData")
load("D:/我的课题/RA/0321上传/RA_CD4T_labeled_0328.RData")
load("D:/我的课题/RA/0321上传/RA_DC_labeled_0201.RData")
# load("C:/Users/aking/Desktop/0321?洗?/pbmc_CD8NK_0214.RData")
#load("C:/Users/aking/Desktop/0321?洗?/RA_DC_and_markers_labeled_0215.RData")


head(RA_DC@meta.data)
head(RA_CD4T@meta.data)
head(Mk@meta.data)
head(B@meta.data)

table(B@meta.data$tissue)
table(Mk@meta.data$tissue)
table(RA_CD4T@meta.data$samples)
table(RA_DC@meta.data$samples)

Idents(B)  %>% table()
Idents(Mk)  %>% table()
Idents(RA_CD4T)  %>% table()
Idents(RA_DC) %>% table()



Idents(B) <- "tissue"
PBMC_B <- subset(B, idents = "PBMC")
Idents(PBMC_B) <- "cell.type"


# Idents(Mk) <- "tissue"
# PBMC_Mk_mono <- subset(Mk, idents = "PBMC")
# Idents(PBMC_Mk_mono) <- "cell.type"

Idents(RA_DC) <- "tissue"
PBMC_DC <- subset(RA_DC, idents="PBMC")
Idents(PBMC_DC) <- "cell.type"
table(PBMC_DC@meta.data$cell.type)


# PBMC_DC %>% Idents() %>% table()
PBMC_B %>% Idents() %>% table()
# PBMC_Mk_mono %>% Idents() %>% table()
RA_CD4T %>% Idents() %>% table()

#Idents(RA_CD4T) <- "cell.type"
#RA_CD4T %>% Idents() %>% table()

merged_seurat <- merge(x = PBMC_B, y = list(RA_CD4T, PBMC_DC))

merged_seurat %>% Idents() %>% table()

table(merged_seurat@meta.data$samples)

Sub_merged__seurat <- subset(merged_seurat, 
                        idents= c("HLA-DRB5+ Memory B",
                                  "HLA-DRB5+ DC", "HLA-DRB5- DC", "Effector CD4T_HLAhiGZMB+",
                                  "Effector CD4T_S100A8hiGZMB+","Effectot CD4T_GZMB+", "Effectot CD4T_GZMK+"))


Sub_merged__seurat
Sub_merged__seurat %>% Idents() %>% table()

##########
dir.create("../Merged_PBMC_data")
setwd("../Merged_PBMC_data")
Idents(Sub_merged__seurat) <- "cell.type"
Sub_merged__seurat

# take raw data and normalise it
count_raw <- GetAssayData(Sub_merged__seurat, slot = "counts", assay = "RNA")
count_raw[1,1:55]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
count_norm[1:5,1:55]
write.table(count_norm, file = "cellphonedb_count.txt", sep="\t", quote=F)

# generating meta file
Cell <- rownames(Sub_merged__seurat@meta.data)
meta_data <- cbind(Cell, Sub_merged__seurat@meta.data[,"cell.type", drop=F])
head(meta_data)
# meta_data$cell_type <- meta_data$cell.type
head(meta_data)
write.table(meta_data, file = "cellphonedb_meta.txt", sep="\t", quote=F, row.names=F)







#####3














Idents(Sub_merged__seurat) <- "ACPA"
head(Sub_merged__seurat@meta.data)
table(Sub_merged__seurat@meta.data$ACPA)

HC_Sub_merged__seurat <- subset(Sub_merged__seurat, idents = "HC")
Idents(HC_Sub_merged__seurat) <- "cell.type"
VlnPlot(HC_Sub_merged__seurat, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 2)
HC_Sub_merged__seurat

Positive_Sub_merged__seurat <- subset(Sub_merged__seurat, idents = "Positive")
Idents(Positive_Sub_merged__seurat) <- "cell.type"
Positive_Sub_merged__seurat
#############################







#############################
#############################
###make metadata and count data file for cellphoneDB analysi
dir.create("./HC_pbmc_data")
setwd("HC_pbmc_data")


HC_Sub_merged__seurat <- subset(HC_Sub_merged__seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 3500)
Idents(HC_Sub_merged__seurat) <- "cell.type"
HC_Sub_merged__seurat


RA_metadata <- HC_Sub_merged__seurat@meta.data[,c(22,23)]
head(RA_metadata)
RA_metadata$Cell <- rownames(RA_metadata) 
RA_metadata$cell_type <- RA_metadata$cell.type

RA_metadata <- RA_metadata[,c(3,4)]
head(RA_metadata)
RA_metadata <- as.data.frame(RA_metadata[,c(1,2)])
head(RA_metadata)
write.table(RA_metadata, file="RA_seurat_metadata.txt", row.names = FALSE, quote = F, sep = "\t")

##RNA expression count file
# RA_seurat_count <- HC_Sub_merged__seurat@assays$RNA[,]
# RA_seurat_count[1:5,1:55]
# write.table(RA_seurat_count, file = "RA_seurat_count.txt",sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
# HC_Sub_merged__seurat

# take raw data and normalise it
count_raw <- GetAssayData(HC_Sub_merged__seurat, slot = "counts", assay = "RNA")
count_raw[1,1:55]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
count_norm[1:5,1:55]
write.table(count_norm, file = "cellphonedb_count.txt", sep="\t", quote=F)

# generating meta file
Cell <- rownames(HC_Sub_merged__seurat@meta.data)
meta_data <- cbind(Cell, HC_Sub_merged__seurat@meta.data[,"cell.type", drop=F])
head(meta_data)
# meta_data$cell_type <- meta_data$cell.type
head(meta_data)
write.table(meta_data, file = "cellphonedb_meta.txt", sep="\t", quote=F, row.names=F)

#############################
#############################
#############################
###make metadata and count data file for cellphoneDB analysi
dir.create("../Nega_PBMC_data")
setwd("../Nega_PBMC_data")
merged_seurat_sub <- subset(Sub_merged__seurat, idents="Negative")
Idents(merged_seurat_sub) <- "cell.type"
merged_seurat_sub

# take raw data and normalise it
count_raw <- GetAssayData(merged_seurat_sub, slot = "counts", assay = "RNA")
count_raw[1,1:55]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
count_norm[1:5,1:55]
write.table(count_norm, file = "cellphonedb_count.txt", sep="\t", quote=F)

# generating meta file
Cell <- rownames(merged_seurat_sub@meta.data)
meta_data <- cbind(Cell, merged_seurat_sub@meta.data[,"cell.type", drop=F])
head(meta_data)
# meta_data$cell_type <- meta_data$cell.type
head(meta_data)
write.table(meta_data, file = "cellphonedb_meta.txt", sep="\t", quote=F, row.names=F)



# RA_metadata <- Negative_Sub_merged__seurat@meta.data[,c(22,23)]
# head(RA_metadata)
# RA_metadata$Cell <- rownames(RA_metadata) 
# RA_metadata$cell_type <- RA_metadata$cell.type
# RA_metadata <- RA_metadata[,c(3,4)]
# head(RA_metadata)
# RA_metadata <- as.data.frame(RA_metadata[,c(1,2)])
# write.table(RA_metadata, file="RA_seurat_metadata.txt", row.names = FALSE, quote = F, sep = "\t")

##RNA expression count file
# RA_seurat_count <- Negative_Sub_merged__seurat@assays$RNA[,]
# write.table(RA_seurat_count, file = "RA_seurat_count.txt",sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
# Negative_Sub_merged__seurat

#############################
#############################
#############################
###make metadata and count data file for cellphoneDB analysi
dir.create("../Pos_PBMC_data")
setwd("../Pos_PBMC_data")
merged_seurat_sub <- subset(Sub_merged__seurat, idents="Negative")
Idents(merged_seurat_sub) <- "cell.type"
merged_seurat_sub

# take raw data and normalise it
count_raw <- GetAssayData(merged_seurat_sub, slot = "counts", assay = "RNA")
count_raw[1,1:55]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
count_norm[1:5,1:55]
write.table(count_norm, file = "cellphonedb_count.txt", sep="\t", quote=F)

# generating meta file
Cell <- rownames(merged_seurat_sub@meta.data)
meta_data <- cbind(Cell, merged_seurat_sub@meta.data[,"cell.type", drop=F])
head(meta_data)
# meta_data$cell_type <- meta_data$cell.type
head(meta_data)
write.table(meta_data, file = "cellphonedb_meta.txt", sep="\t", quote=F, row.names=F)









# RA_metadata <- Positive_Sub_merged__seurat@meta.data[,c(22,23)]
# head(RA_metadata)
# RA_metadata$Cell <- rownames(RA_metadata) 
# RA_metadata$cell_type <- RA_metadata$cell.type
# RA_metadata <- RA_metadata[,c(3,4)]
# head(RA_metadata)
# RA_metadata <- as.data.frame(RA_metadata[,c(1,2)])
# write.table(RA_metadata, file="RA_seurat_metadata.txt", row.names = FALSE, quote = F, sep = "\t")

##RNA expression count file
# RA_seurat_count <- Positive_Sub_merged__seurat@assays$RNA[,]
# write.table(RA_seurat_count, file = "RA_seurat_count.txt",sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
# Positive_Sub_merged__seurat


####################test demo---it worked
merged_seurat <- merge(x = Mk, y = list(B, RA_CD4T))

MK_1 <- subset(Mk, cells = 1:10)
B_1 <- subset(B,  cells=1:10)
RA_CD4T_1 <- subset(RA_CD4T, cells = 1:10)

merged_seurat <- merge(x = MK_1, y = list(B_1, RA_CD4T_1))


head(merged_seurat@meta.data)
#######
RA_metadata <- merged_seurat@meta.data[,c(22,23)]
head(RA_metadata)
RA_metadata$Cell <- rownames(RA_metadata) 
RA_metadata$cell_type <- RA_metadata$cell.type

RA_metadata <- RA_metadata[,c(3,4)]
head(RA_metadata)
RA_metadata <- as.data.frame(RA_metadata[,c(1,2)])
write.table(RA_metadata, file="RA_seurat_metadata.txt", row.names = FALSE, quote = F, sep = "\t")


#####
RA_seurat_count <- merged_seurat@assays$RNA[,]

#RA_seurat_count[,1:5]
class(RA_seurat_count)

write.table(RA_seurat_count, file = "RA_seurat_count.txt",sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)