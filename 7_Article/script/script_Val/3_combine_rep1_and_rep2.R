# load the Seurat object
path_RNA_combined  <- "7_Article/results/tab_seurat_correct/RNA_combined.rds"
RNA_combined <- readRDS(path_RNA_combined)

#path_RNA_rep1  <- "7_Article/results/tab_seurat_correct/RNA_rep1.rds"
#RNA_rep1 <- readRDS(path_RNA_rep1)

#path_RNA_rep2  <- "7_Article/results/tab_seurat_correct/RNA_rep2.rds"
#RNA_rep2 <- readRDS(path_RNA_rep2)


# Joindre les couches pour obtenir une seule matrice d'expression
combined_seurat_join[["RNA"]] <- JoinLayers(combined_seurat)

# Sauvegarder l'objet Seurat
saveRDS(combined_seurat_join, "7_Article/results/tab_seurat_correct/RNA_combined_join.rds")

