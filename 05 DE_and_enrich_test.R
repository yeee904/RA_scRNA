## 
### ---------------
###
### Create: Yi Liu
### Date: 2020/08/03
### Update Log: 2020/08/04  
### This script aims to find differential genes and do enrich test;calculate module score.
### only an example
### ---------------

#DE gene test using  Seurat "FindMarkers" function
library(Seurat)
GNLY.memsw.markers <- FindMarkers(B,ident.1 = "GNLY+ Mem-sw B",ident.2 = "Mem-sw B",logfc.threshold = 0)
B.markers <- FindAllMarkers(B,logfc.threshold = 0.5,only.pos = T)

###enrichment analysis using clusterprofiler
library(clusterProfiler)
library(KEGG.db)
library(dplyr)
my_enrichKEGG <- function(de_clusters,...){
  my_kegg <- function(de_clusters){
    genes <- de_clusters$entrezid
    kegg_result <- enrichKEGG(genes,...)
    return(kegg_result)
  }
  #data_kegg <- plyr::dlply(.data=de_clusters,.variable=.(cluster),my_kegg)
  pieces <- split(de_clusters,list(de_clusters$cluster))
  data_kegg <- lapply(pieces, my_kegg)
  get_result <- function(x){
    result=x@result
    return(result)
  }
  combine_result <- plyr::ldply(.data = data_kegg,.fun =get_result)
  combine_result$cluster <- combine_result$`.id`
  return(combine_result)
}
my_enrichGO <- function(de_clusters,...){
  my_go <- function(de_clusters){
    genes <- de_clusters$entrezid
    go_result <- enrichGO(genes,...)
    return(go_result)
  }
  #data_kegg <- plyr::dlply(.data=de_clusters,.variable=.(cluster),my_kegg)
  pieces <- split(de_clusters,list(de_clusters$cluster))
  data_go <- lapply(pieces, my_go)
  get_result <- function(x){
    result=x@result
    return(result)
  }
  combine_result <- plyr::ldply(.data = data_go,.fun =get_result)
  combine_result$cluster <- combine_result$`.id`
  return(combine_result)
}
macro.markers <- FindAllMarkers(macro,logfc.threshold = 0.5,only.pos = T)
macro.markers <- subset(macro.markers, p_val_adj<0.05)
macro.markers %>% filter(!(grepl("^RP",macro.markers$gene) | grepl("^MT-",macro.markers$gene)))->macro.markers
entrez_genes <- bitr(macro.markers$gene, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Hs.eg.db")
head(entrez_genes)
macro_de_clusters <- subset(macro.markers,macro.markers$gene %in% entrez_genes$SYMBOL)
macro_de_clusters$entrezid <- plyr::mapvalues(x=macro_de_clusters$gene,from = entrez_genes$SYMBOL,to=entrez_genes$ENTREZID)
head(macro_de_clusters)
macro_kegg <- my_enrichKEGG(macro_de_clusters,
                            organism='hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            use_internal_data=T)
macro_GO <- my_enrichGO(macro_de_clusters,
                        OrgDb="org.Hs.eg.db",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
##calculate module score
genes=subset(genes,genes %in% rownames(RA_DC))
genes=list(genes)
seurat_obj <- AddModuleScore(object = seurat_obj,features = genes,name = score)