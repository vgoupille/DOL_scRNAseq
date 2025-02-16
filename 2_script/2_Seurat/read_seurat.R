# -*- coding: utf-8 -*-


#singularity exec /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/6_Containers/seurat.sif Rscript 2_script/2_Seurat/read_seurat.R

# Charger les bibliothèques nécessaires
library(Seurat)
library (ggplot2)
library(dplyr)
library(patchwork)
# Read Seurat object

#4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds
#4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds

#DOL_scRNAseq_filtered  <- readRDS("4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds")

#DOL_scRNAseq_raw  <- readRDS("4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds")
pmbc <- readRDS("4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds")

# Print head of the object
#head(DOL_scRNAseq_filtered, 100)
#head(DOL_scRNAseq_raw, 100) # je veux les 100 premières lignes
head(pmbc, 100) # je veux les 100 premières lignes

# Attaching SeuratObject
#                            orig.ident nCount_RNA nFeature_RNA
# AAACATCG_AAACATCG_AAACGATA   AAACATCG          3            3
# AACGTGAT_AAACATCG_AAACGATA   AACGTGAT          4            3
# AAACATCG_AACGTGAT_AAACGATA   AAACATCG          6            6
# AACGTGAT_AACGTGAT_AAACGATA   AACGTGAT          1            1
# AAACATCG_AAACATCG_ACTCGTAA   AAACATCG        150           21
# AACGTGAT_AAACATCG_ACTCGTAA   AACGTGAT        260           20
# AAACATCG_AACGTGAT_ACTCGTAA   AAACATCG        480           36
# AACGTGAT_AACGTGAT_ACTCGTAA   AACGTGAT        552           32
# AAACATCG_AAACATCG_CATGATCA   AAACATCG        763           60
# AACGTGAT_AAACATCG_CATGATCA   AACGTGAT       1353           66
# AAACATCG_AACGTGAT_CATGATCA   AAACATCG       1081           65
# AACGTGAT_AACGTGAT_CATGATCA   AACGTGAT       1023           54
# AAACATCG_AAACATCG_CTGCTTTG   AAACATCG       7535          348
# AACGTGAT_AAACATCG_CTGCTTTG   AACGTGAT       9007          349
# AAACATCG_AACGTGAT_CTGCTTTG   AAACATCG      23664          394
# AACGTGAT_AACGTGAT_CTGCTTTG   AACGTGAT      27104          398

# il y a 3 colonnes dans le tableau, orig.ident, nCount_RNA, nFeature_RNA
# orig.ident est le nom de la cellule
# nCount_RNA est le nombre de lecture pour cette cellule
# nFeature_RNA est le nombre de gène pour cette cellule

# dans cette experimentation seulement 16 "cellules" => regroupe plusieurs cellules à cause de l'histoire des barcodes

# Print the object info 
pmbc
#An object of class Seurat 
#6248 features across 16 samples within 1 assay 
#Active assay: RNA (6248 features, 0 variable features)

# Lets examine a few genes in the first thirty cells
#pmbc[c ("NfuA", "soxR"), 1:16]