## 
### ---------------
###
### Create: Shanzhao Jin
### Date: 2020/04/25
### Update Log: 2020/08/04  
### This script aims to predict cell-cell interaction using cellphoneDB.
### 
### ---------------


#!/bin/usr/sh
source /home/jinshanzhao/cpdb-venv/bin/activate
# RA_seurat_metadata.txt = "$1"
# RA_seurat_count.tx = "$2"

mkdir mycellphoneda
#Example with running the statistical method
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name --threads=5 --output-path=mycellphoneda --threshold=0.1

#Example without using the statistical method
# cellphonedb method analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name --output-path=mycellphoneda

#Plotting statistical method results
cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --output-path=./mycellphoneda 
cellphonedb plot heatmap_plot cellphonedb_meta.txt --pvalues-path=./mycellphoneda/pvalues.txt --output-path=./mycellphoneda
# cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows mycellphoneda/rows.txt --columns mycellphoneda/columns.txt --output-path=./mycellphoneda --output-name=selected_dotplot.pdf

#Subsampling
# cellphonedb method analysis cellphonedb_meta.txt cellphonedb_count.txt --subsampling --subsampling-log false --subsampling-num-cells 3000


cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/rows1.txt --columns /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/columns1.txt --output-path=./mycellphoneda --output-name=selected_dotplot_1.pdf

cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/rows2.txt --columns /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/columns2.txt --output-path=./mycellphoneda --output-name=selected_dotplot_2.pdf

cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/rows3.txt --columns /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/columns3.txt --output-path=./mycellphoneda --output-name=selected_dotplot_3.pdf

cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/rows4.txt --columns /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/columns4.txt --output-path=./mycellphoneda --output-name=selected_dotplot_4.pdf

cellphonedb plot dot_plot --means-path=./mycellphoneda/means.txt --pvalues-path=./mycellphoneda/pvalues.txt --rows /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/rows5.txt --columns /Bailab7/PROJECT/jinshanzhao/singlecell/RA/result/cellphoneDB/RA_result/RA_cellphoneDB-20200507/PBMC-cellDB/Selected_pairs/columns5.txt --output-path=./mycellphoneda --output-name=selected_dotplot_5.pdf


#Subsampling
