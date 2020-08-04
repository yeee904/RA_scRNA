# This script is designed to pipeline Seurat object output into Monocle3
# This script also contains all the necessary functionalities in Monocle3 as 
# This could also convert a 2D Seurat object and visualize/analyze it in 3D
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# This script is originally written for local machines but adaptations have also been included in annotations


### Require:: 'BiocManager', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'reticulate', 'htmlwidgets'

# This is a required python package
#reticulate::py_install("louvain")

# This is installing the actual monocle3
#devtools::install_github('cole-trapnell-lab/monocle3')

cell_type_colors <- c("#A6CEE3", "#1F78B4",
                      "#fadc58","#e2e589",
                      "#B2DF8A", "#33A02C",
                      "#f47d42","#FF6347",
                      "#CAB2D6", "#6A3D9A","#8d798a", "#5d4e60","#74a57f")



### Installing the packages

library(Seurat)
library(monocle3)
library(htmlwidgets)


### Specifying input parameters

# args = commandArgs(trailingOnly = T)
# input.dir <- args[1]
# output.dir.input <- args[2]
# upn <- args[3]
# Dim <- args[4]
# Celltypefile <- args[5]
# nPC <- args[6]
# cluster.res <- args[7]
# output.dir <- sprintf("%s/%s.%s.Seurat.to.Monocle3.%s", output.dir.input, upn, date, Dim)

date = gsub("2020-","2020",Sys.Date(),perl=TRUE);
date = gsub("-","",date)
Dim = "2D"
upn = "B_cell"
input.dir <- getwd()
nPC <- 30
cluster.res <- 0.5
output.dir <- sprintf("%s/%s.%s.Seurat.to.Monocle3.%s", getwd(), upn, date, Dim)
dir.create(file.path(output.dir), showWarnings = FALSE)


### Reading in Seurat object

print("Readingin Seurat objects")

#seurat <- readRDS(file = input.dir)
# load("B_subset")
Idents(B) <- "tissue"

B_PBMC <- subset(B, idents="PBMC")
Idents(B_PBMC) <- "cell.type"

seurat <- B_PBMC
DefaultAssay(seurat) <- "integrated"
### Re-dimension reduction for 3D rendering
if (Dim = "3D"){
  print ("Running UMAP 3D")
  
  seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:nPC, n.components = 2)
  print("Clustering 3D")
  seurat <- FindNeighbors(object=seurat, dims=1:nPC)
 # seurat <- FindClusters(object=seurat, resolution=cluster.res)
  seurat[[sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)]] <- Idents(object = seurat)
}

DimPlot(seurat)
DimPlot(seurat, cols = cell_type_colors)
### Building the necessary parts for a basic cds

# part 1, gene annotations

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), 
                                 row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part 2, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
head(cell_metadata)
# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
expression_matrix[1:5,1:5]

### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


### Could be a space-holder, but essentially fills out louvain parameters

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds_from_seurat@reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings


### Assign feature loading for downstream module analysis

cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings


### Learn graph, this step usually takes a significant period of time for larger samples

print("Learning graph, which can take a while depends on the sample")

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)


### Reading in cell type information

print("Reading in cell type information")

celltype_infor <- B_PBMC@meta.data[,c(7,21,22,23)]
celltype_infor$cell.barcodes.i. <- rownames(celltype_infor)
# celltype_infor <- celltype_infor[,2:3]
head(celltype_infor)

# write(celltype_infor, file = "celltype_infor_predictions.csv")
# celltype_infor <- "celltype_infor_predictions.csv"


#lineage.table <- read.table(file = celltype_infor, sep = "", header = T, row.names = F)
#lineage.table$cell.barcodes.i. <- paste(lineage.table$cell.barcodes.i., "-1", sep = "")
#lineage.table <- lineage.table[, 1:2]
#lineage.table$lineage1 <- gsub("_\\d+","", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("TCELLA","T-CELL", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("BCELLA","B-CELL", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("NKA","NK", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("DENDRA","DENDRITIC", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("DENDA","DENDRITIC", lineage.table$lineage1)
#lineage.table$lineage1 <- gsub("[0-9]+", "", lineage.table$lineage1)
#lineage.table$cell.barcodes.i. <- gsub("-1", "", lineage.table$cell.barcodes.i.)

#indices <- match(pData(cds_from_seurat)$barcode, lineage.table$cell.barcodes.i.)
#lineage.table.refined <- lineage.table[indices,]
#lineage.table.final <- lineage.table.refined[, 2]
#colData(cds_from_seurat)$celltype <- lineage.table.final


colData(cds_from_seurat)$celltype <- celltype_infor$cell.type
colData(cds_from_seurat)$ACPA <- celltype_infor$ACPA
colData(cds_from_seurat)$samples <- celltype_infor$samples
colData(cds_from_seurat)$tissue <- celltype_infor$tissue



### Plot cluster info with trajectory

print("Plotting clusters")

if (Dim == "2D") {
  pdf(sprintf("%s/clusters.with.trajectory.%s.pdf", output.dir, Dim), width = 6, height = 6)
  clus <- plot_cells(cds_from_seurat, 
                     color_cells_by = 'cluster',
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE)+
    scale_color_manual(values=cell_type_colors)+
    theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = rep(1))) + 
    theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold",colour = "blue")) 
  print(clus)
  dev.off()
}

if (Dim == "3D") {
  cluster <- plot_cells_3d(cds_from_seurat, 
                           color_cells_by = 'cluster') +
    scale_color_manual(values=cell_type_colors)+
    theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = rep(1))) + 
    theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  saveWidget(cluster, file = sprintf("%s/clusters.with.trajectory.%s.html", output.dir, Dim))
}

if (Dim == "3D") {
  cluster <- plot_cells_3d(cds_from_seurat, 
                           color_cells_by = 'celltype') +
    scale_color_manual(values=cell_type_colors) +
    theme(axis.text = element_text(size = 10,face = "bold"), 
          axis.title = element_text(size = rep(1))) + 
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  htmlwidgets::saveWidget(cluster, file = sprintf("%s/celltype.with.trajectory.%s.html", output.dir, Dim))
}

### Plot cell type info with trajectory

print("Plotting cell type info")

if (Dim == "2D") {
  pdf(sprintf("%s/celltype.with.trajectory.%s.pdf", output.dir, Dim), width = 6, height = 5.5)
  ctype <- plot_cells(cds_from_seurat, 
                      color_cells_by = 'celltype',
                      label_groups_by_cluster= FALSE,
                    #  label_cell_groups = FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE,
                      group_label_size = 4) +
    scale_color_manual(values=cell_type_colors)+
    theme(axis.text = element_text(size = 10,face = "bold"), 
          axis.title = element_text(size = rep(10))) + 
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  print(ctype)
  dev.off()
}

if (Dim == "2D") {
  pdf(sprintf("%s/ACPA.with.trajectory.%s.pdf", output.dir, Dim), width = 6, height = 5)
  ctype <- plot_cells(cds_from_seurat, 
                      group_label_size = 3,
                      color_cells_by =  "ACPA",
                      label_cell_groups = FALSE,
                      label_groups_by_cluster=FALSE,
                      label_branch_points=FALSE,
                      label_leaves=FALSE) +
    scale_color_manual(values=cell_type_colors) +
    theme(axis.text = element_text(size = 10,face = "bold"),
          axis.title = element_text(size = rep(10))) + 
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  print(ctype)
  dev.off()
}

if (Dim == "2D") {
  pdf(sprintf("%s/tissue.with.trajectory.%s.pdf", output.dir, Dim), width = 6, height = 5)
  ctype <- plot_cells(cds_from_seurat, 
                      group_label_size = 3,
                      color_cells_by =  "tissue",
                      label_cell_groups = FALSE,
                      label_groups_by_cluster=FALSE,
                      label_branch_points=FALSE,
                      label_leaves=FALSE) +
    scale_color_manual(values=cell_type_colors) +
    theme(axis.text = element_text(size = 10,face = "bold"),
          axis.title = element_text(size = rep(10))) + 
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  print(ctype)
  dev.off()
}

if (Dim == "3D") {
  celltype <- plot_cells_3d(cds_from_seurat, 
                            color_cells_by = 'celltype',
                            label_groups_by_cluster=FALSE,
                            label_leaves=FALSE,
                            label_branch_points=FALSE) +
    scale_color_manual(values=cell_type_colors)+
    theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = rep(3))) + 
    theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold",colour = "blue"))
  saveWidget(celltype, file = sprintf("%s/celltype.with.trajectory.%s.html", output.dir, Dim))
}


### You can also plot out different subclone information, which is not shown here


# Before we order the cells in pseudotime, we need to identify which is the pricipal node/cells
# We can either use their default GUI (which unfortunately is not compatible with 3D), or in the case of NBM, select the HSCs
# For the sake of completeness, we will feed in the cells they are labeled as HSCs, but in other samples for example HSCs, it could be a different story
# plus I am not sure if the GUI will show up on the cluster

# This is using the GUI
# cds_from_seurat <- order_cells(cds_from_seurat)

root_cell_list_1 <- grep(c("Naive B"), colData(cds_from_seurat)$celltype)
#root_cell_list_2 <- grep("CCL+ memory B" , colData(cds_from_seurat)$celltype)

#root_cell_list_2

root_cell_list <- counts(cds_from_seurat)@Dimnames[[2]][root_cell_list_1]

## Defined a function which calculate a "celltype" as start point in pseudotime
get_earliest_principal_node <- function(cds_from_seurat, time_bin=c("Naive B")){
  cell_ids <- which(colData(cds_from_seurat)[, "celltype"] == time_bin)
  closest_vertex <- cds_from_seurat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_from_seurat), ])
  root_pr_nodes <-igraph::V(principal_graph(cds_from_seurat)[["UMAP"]])$name[as.numeric(names
                                                                            (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

cds_from_seurat = order_cells(cds_from_seurat, root_pr_nodes=get_earliest_principal_node(cds_from_seurat))



# Sometimes, in this case, the root cells are dispered so single cell has to be identified to input into the order cells function
# cds_from_seurat <- order_cells(cds_from_seurat, root_cells = "CTTGATTCACAAAGTA")
# cds_from_seurat <- order_cells(cds_from_seurat, root_cells = root_cell_list)


### Now we plot pseudotime

print("Plotting pseudotime")

plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) +  
  theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = rel(1))) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold",colour = "blue")) +
  ggsave("B_subset_graph_by_auto_ordered_Naive_B.png", width=6, height=5, dpi = 600)

plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) +  
  facet_grid(~ACPA)+
  theme(axis.text = element_text(size = 10,face = "bold"), axis.title = element_text(size = rel(1))) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold",colour = "blue")) +
  ggsave("B_subset_graph_by_auto_ordered_Naive_B_facet.png", width=10, height=4, dpi = 600)

if (Dim == "2D") {
  pdf(sprintf("%s/pseudotime.with.trajectory.%s.pdf", output.dir, Dim), width = 6, height = 5)
  ptime <- plot_cells(cds_from_seurat, 
                      color_cells_by = 'pseudotime',
                      label_cell_groups = FALSE,
                      label_groups_by_cluster= TRUE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE)
  print(ptime)
  dev.off()
}

if (Dim == "3D") {
  pseudotime <- plot_cells_3d(cds_from_seurat, 
                              color_cells_by = 'pseudotime',
                              label_groups_by_cluster=FALSE,
                              label_leaves=FALSE,
                              label_branch_points=FALSE)
  saveWidget(pseudotime, file = sprintf("%s/pseudotime.with.trajectory.%s.html", output.dir, Dim))
}


### Now we can subset the cells and do branch expression/pseudotime analysis
### This first part is to identify genes that are differentially expressed on the different paths through the trajectory

# This function is disabled for now

p <- 0
if (p > 0) {
  
  print("calculating most significant genes across the whole trajectory, EXRREMELY time-consuming")
  
  cds_from_seurat_pr_test_res = graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=4)
  pr_deg_ids = row.names(subset(cds_from_seurat_pr_test_res, q_value < 0.05))
  
  gene_list_1 <- pr_deg_ids[1:5]
  
  
  ### Plotting individual genes on the map
  
  print("Plotting individual genes across the whole trajectory")
  
  gene_list_1 = c("IRF8","GPR183","MZB1","CCL18","GNLY","HLA-DRB5","IGHM","IGHA1","IGHG1","IGHG3")
  
  
  if (Dim == "2D") {
    pdf(sprintf("%s/top4.highly.sig.genes.with.trajectory.%s.pdf", output.dir, Dim), width = 8, height = 6)
    hgenes <- plot_cells(cds_from_seurat, group_label_size = 10,
                         genes = gene_list_1,
                         show_trajectory_graph=FALSE,
                         label_cell_groups=FALSE,
                         label_leaves=FALSE) +
      theme(axis.text = element_text(size = 10,face = "bold"), 
            axis.title = element_text(size = rep(10))) + 
      theme(legend.text = element_text(size = 10), 
            legend.title = element_text(size = 10, face = "bold",colour = "blue"))
    print(hgenes)
    dev.off()
  }
  
  if (Dim == "3D") {
    indi_genes <- plot_cells_3d(cds_from_seurat, 
                                genes = gene_list_1,
                                label_groups_by_cluster=FALSE,
                                label_leaves=FALSE,
                                label_branch_points=FALSE)
    saveWidget(indi_genes, file = sprintf("%s/top4.highly.sig.genes.with.trajectory.%s.html", output.dir, Dim))
}
  
  ### Now we can collect the trajectory-variable genes into modules for co-reg elements
  My_genes_cds_pr_test_res = graph_test(cds_from_seurat,  neighbor_graph="principal_graph", cores=1)
  pr_deg_ids = row.names(subset(My_genes_cds_pr_test_res, q_value < 0.05))
  gene_module_df = monocle3:::find_gene_modules(cds_from_seurat[pr_deg_ids,])
  # gene_module_df = monocle3:::find_gene_modules(cds_from_seurat[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
  
  cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat)),
                                 cell_group=colData(cds_from_seurat)$celltype)
  agg_mat = aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
  row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
  pheatmap::pheatmap(agg_mat,
                     scale="column", clustering_method="ward.D2")
  pheatmap::pheatmap(agg_mat,
                     scale="column", clustering_method="ward.D2",
                     fontsize=6, width=5, height=8,
                     file="B_suset_pseudotime_module_heatmap.png")
  
  
  ### Now we can plot these modules
  
  if (Dim == "2D") {
    pdf(sprintf("%s/sig.modules.with.trajectory.%s.pdf", output.dir, Dim), width = 14, height = 10)
    smodules <- plot_cells(cds_from_seurat,
                           genes=gene_module_df %>% filter(module %in% c(7, 18, 21, 10,17,23,6,15,12,9)),
                           label_cell_groups=FALSE, 
                           show_trajectory_graph=FALSE) +
      theme(axis.text = element_text(size = 10,face = "bold"), 
            axis.title = element_text(size = rep(10))) + 
      theme(legend.text = element_text(size = 10), 
            legend.title = element_text(size = 10, face = "bold",colour = "blue"))
    print(smodules)
    dev.off()
  }
  
  if (Dim == "3D") {
    indi_modules <- plot_cells_3d(cds_from_seurat, 
                                  genes=gene_module_df %>% filter(module %in% c(7)))
    htmlwidgets::saveWidget(indi_modules, file = sprintf("%s/sig.modules.with.trajectory_module7.%s.html", output.dir, Dim))
  }
 
indi_modules <- plot_cells_3d(cds_from_seurat, 
                                genes=gene_module_df %>% filter(module %in% c(18)))
htmlwidgets::saveWidget(indi_modules, 
                        file = sprintf("%s/sig.modules.with.trajectory_module18.%s.html", 
                                       output.dir, Dim))
indi_modules <- plot_cells_3d(cds_from_seurat, 
                              genes=gene_module_df %>% filter(module %in% c(21)))
htmlwidgets::saveWidget(indi_modules, 
                        file = sprintf("%s/sig.modules.with.trajectory_module21.%s.html", 
                                       output.dir, Dim))  
  
   
  ### Or we can plot them against pseudotime
  # Say you want to study the expression of certain gene (MPO) in certain lineages (HSC - MEP - CEP - Monocyte)
  
  Naive_genes = c("IGHM")
  Mono_lineage_cds_from_seurat = cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% Naive_genes,
                                                 colData(cds_from_seurat)$celltype %in% c("Naive B")]
  Mono_lineage_cds_from_seurat
  print("Plotting genes against pseudotime")
  
  pdf(sprintf("%s/Genes.against.pseudotime.%s.pdf", output.dir, Dim), width = 20, height = 20)
  plot_genes_in_pseudotime(Mono_lineage_cds_from_seurat,
                           color_cells_by="celltype",
                           min_expr=0.5)
  dev.off()
  
}


### Now we can subset the cells and do branch expression/pseudotime analysis
### This second part is to identify genes that are differentially expressed only within this selected region but not necessarily with the trajectories

# This function is disabled for now as well

if (p > 1) {
  
  # Subsetting the cells
  # Yet this step required GUI, which I am not sure how it will pan out on the cluster
  
  cds_from_seurat_subset = choose_cells(cds_from_seurat)
  
  
  # Identify the significant genes
  
  print("Subset analysis")
  
  subset_pr_test_res = graph_test(cds_from_seurat_subset, cores=4)
  pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < 0.05))
  
  
  # Group them into modules
  # Here the resolution parameters and k might vary so feel free to adjust them as well
  # Sometimes when there are too few genes we might want to 
  
  gene_module_df = monocle3:::find_gene_modules(cds_from_seurat_subset[pr_deg_ids,], resolution=1e-5, verbose = T, k = 10)
  
  
  # Plot the scores for the modules
  cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat_subset)), cell_group=clusters(cds_from_seurat_subset)[colnames(cds_from_seurat_subset)])
  agg_mat = aggregate_gene_expression(cds_from_seurat_subset, gene_module_df, cell_group_df)
  module_dendro = hclust(dist(agg_mat))
  gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])
  
  
  # Plot the modules along the trajectory
  if (Dim == "2D") {
    pdf(sprintf("%s/subset.sig.modules.with.trajectory.%s.pdf", output.dir, Dim), width = 12, height = 10)
    submodules <- plot_cells(cds_from_seurat,
                             genes=gene_module_df,
                             label_cell_groups=FALSE,
                             show_trajectory_graph=FALSE)
    scale_color_manual(values=cell_type_colors)+
      theme(axis.text = element_text(size = 10,face = "bold"), 
            axis.title = element_text(size = rep(10))) + 
      theme(legend.text = element_text(size = 10), 
            legend.title = element_text(size = 10, face = "bold",colour = "blue"))
    print(submodules)
    dev.off()
  }
  
  if (Dim == "3D") {
    subset_indi_modules <- plot_cells_3d(cds_from_seurat, 
                                         genes=gene_module_df,
                                         label_groups_by_cluster=FALSE,
                                         label_leaves=FALSE,
                                         label_branch_points=FALSE)
    saveWidget(subset_indi_modules, file = sprintf("%s/subset.sig.modules.with.trajectory.%s.html", output.dir, Dim))
  }
}

### End of the script

