

# Import the libraries
library (Seurat)

# load the Seurat object

path_RNA_combined  <- "7_Article/results/tab_seurat_correct/RNA_combined.rds"
RNA_combined <- readRDS(path_RNA_combined)


#method to combine the layers of the Seurat object
# https://satijalab.github.io/seurat-object/reference/SplitLayers.html

Assays(RNA_combined)  # Pour voir les assays disponibles
Layers(RNA_combined) # check the layers


# Joindre les couches pour obtenir une seule matrice d'expression
RNA_combined[["RNA"]] <- JoinLayers(RNA_combined[["RNA"]])


RNA_combined


# Chemin du fichier a sauvegarder
file_path <- "7_Article/results/tab_seurat_correct/RNA_combined_join.rds"



# Sauvegarder l'objet Seurat
saveRDS(RNA_combined, file_path)



# Vérifier si le fichier existe
if (file.exists(file_path)) {
    # Obtenir les informations sur le fichier
    info_fichier <- file.info(file_path)
    
    # Charger l'objet Seurat
    seurat_object_raw <- readRDS(file_path)
    
    # Afficher les informations du fichier
    cat("✓ Fichier créé avec succès :", file_path, "\n")
    cat("Taille du fichier :", round(info_fichier$size/1024/1024, 2), "MB\n")
    cat("Date de création :", format(info_fichier$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Afficher les dimensions de l'objet Seurat
    cat("\nInformations sur l'objet Seurat :\n")
    cat("Nombre de cellules :", ncol(seurat_object_raw), "\n")
    cat("Nombre de gènes :", nrow(seurat_object_raw), "\n")
} else {
    cat("❌ Erreur : Le fichier n'a pas été créé.\n")
}