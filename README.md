# RA_scRNA
Rheumatoid Arthritis project. single cell RNA sequencing data.

Softwares:
1. CellRanger v.2.2.0	
source: 10XGenomics, https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
platform: Linux system(8-core Intel or AMD processor (16 cores recommended); 64GB RAM (128GB recommended); 1TB free disk space;64-bit CentOS/RedHat 6.0 or Ubuntu 12.04)
typical install time: 12h
2. Seurat v3.1.1	(Butler et al., 2018; Stuart et al., 2019).
source: https://satijalab.org/seurat/v3.1/
platform: R, install by "install.packages(Seurat)"
3. clusterProfiler v3.12.0	(Yu et al., 2012)
source: https://guangchuangyu.github.io/software/clusterProfiler/
platform: R, install by "BiocManager::install(clusterPfrofiler)"
4. DAVID	(Huang et al., 2009a, b)
source: https://david.ncifcrf.gov/,online tool,following official structions.
5. Monocle3	(Cao et al., 2019; Qiu et al., 2017a; Qiu et al., 2017b; Trapnell et al., 2014).
source: https://cole-trapnell-lab.github.io/monocle3/
platform:R, install by 
6. CellPhoneDB v2.1.2 (Efremova et al., 2020; Vento-Tormo et al., 2018)
source: https://www.cellphonedb.org/
platform: Python, install by "pip install CellPhoneDB"

Scripts:
01 cellranger.sh
single cell RNA-seq data alignment using CellRanger
output: expression matrix of each sample

02 merge-individual sample.R
merge all individual Seurat objects from the same tissue (PBMC or synovial tissue)

03 extract-major celltype.R
extract major cell types from PBMC and SM.

04 merge-major celltype.R
integrate major cell types from PBMC and SM for further sub-clustering

05 DE_and_enrich_test.R
find differential genes and do enrich test;calculate module score.

06 Seurat-Monocle3.R
using Monocle3 to do trajectory analysis

07-1 samples_prepaired-cellphoneDB.R 07-2 cellphoneDB.sh
predict cell-cell interaction using cellphoneDB. 

