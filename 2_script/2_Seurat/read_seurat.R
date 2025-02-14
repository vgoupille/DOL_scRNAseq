# -*- coding: utf-8 -*-

# Notre repertoire de travail/ Faire un WORKINGDIR
#WORKINGDIR="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/"
#setwd(WORKINGDIR)
#getwd()
 

#Important files which come from the STARsolo output are the ‘barcodes.tsv’, ‘features.tsv’ and ‘UniqueAndMult-Uniform.mtx’ files. These files are used to create the Seurat object. The ‘UniqueAndMult-Uniform.mtx’ file contains the gene expression data in matrix format. The ‘barcodes.tsv’ file contains the cell barcodes and the ‘features.tsv’ file contains the gene names.: 

#Refer to three additional output files for further downstream analysis: the barcode.tsv, features.tsv and UniqueAndMult-Uniform.mtx. The first two files are generated in the ‘Solo.out/GeneFull/’ folder. The UniqueAndMult-Uniform.mtx is located in the ‘Solo.out/ GeneFull/Raw/’ folder.



# Charger les bibliothèques nécessaires
library(Seurat)

# Read Seurat object

#4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds
#4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds

#DOL_scRNAseq_filtered  <- readRDS("4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds")

DOL_scRNAseq_raw  <- readRDS("4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds")


# Print head of the object
#head(DOL_scRNAseq_filtered, 100)
head(DOL_scRNAseq_raw, 100) # je veux les 100 premières lignes

Attaching SeuratObject
                           orig.ident nCount_RNA nFeature_RNA
AAACATCG_AAACATCG_AAACGATA   AAACATCG          3            3
AACGTGAT_AAACATCG_AAACGATA   AACGTGAT          4            3
AAACATCG_AACGTGAT_AAACGATA   AAACATCG          6            6
AACGTGAT_AACGTGAT_AAACGATA   AACGTGAT          1            1
AAACATCG_AAACATCG_ACTCGTAA   AAACATCG        150           21
AACGTGAT_AAACATCG_ACTCGTAA   AACGTGAT        260           20
AAACATCG_AACGTGAT_ACTCGTAA   AAACATCG        480           36
AACGTGAT_AACGTGAT_ACTCGTAA   AACGTGAT        552           32
AAACATCG_AAACATCG_CATGATCA   AAACATCG        763           60
AACGTGAT_AAACATCG_CATGATCA   AACGTGAT       1353           66
AAACATCG_AACGTGAT_CATGATCA   AAACATCG       1081           65
AACGTGAT_AACGTGAT_CATGATCA   AACGTGAT       1023           54
AAACATCG_AAACATCG_CTGCTTTG   AAACATCG       7535          348
AACGTGAT_AAACATCG_CTGCTTTG   AACGTGAT       9007          349
AAACATCG_AACGTGAT_CTGCTTTG   AAACATCG      23664          394
AACGTGAT_AACGTGAT_CTGCTTTG   AACGTGAT      27104          398