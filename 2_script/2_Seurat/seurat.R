# -*- coding: utf-8 -*-

# Notre repertoire de travail/ Faire un WORKINGDIR
#WORKINGDIR="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/"
#setwd(WORKINGDIR)
#getwd()
 

#Important files which come from the STARsolo output are the ‘barcodes.tsv’, ‘features.tsv’ and ‘UniqueAndMult-Uniform.mtx’ files. These files are used to create the Seurat object. The ‘UniqueAndMult-Uniform.mtx’ file contains the gene expression data in matrix format. The ‘barcodes.tsv’ file contains the cell barcodes and the ‘features.tsv’ file contains the gene names.: 

#Refer to three additional output files for further downstream analysis: the barcode.tsv, features.tsv and UniqueAndMult-Uniform.mtx. The first two files are generated in the ‘Solo.out/GeneFull/’ folder. The UniqueAndMult-Uniform.mtx is located in the ‘Solo.out/ GeneFull/Raw/’ folder.

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/filtered/barcodes.tsv

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/filtered/features.tsv

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/raw/UniqueAndMult-Uniform.mtx


#▲ CRITICAL STEP Do not use the ‘matrix.mtx’ folder inside of the ‘Solo.out/GeneFull/’ folder. This file contains only the data from the uniquely aligned genes, whereas we use the data from both unique and multiple alignments.

# Charger les bibliothèques nécessaires
library(Seurat)

#print les fonctions de la bibliothèque Seurat
ls("package:Seurat")

# Définir les chemins des fichiers
barcodes <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/barcodes.tsv"
features <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/features.tsv"
matrix_file <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/raw/UniqueAndMult-Uniform.mtx"

library(Matrix)



print (barcodes)
print (features)
print (matrix_file)


# https://satijalab.org/seurat/reference/#preprocessing

#ReadMtx() : Load in data from remote or local mtx files
# Lire les données
#expression_matrix <- ReadMtx(mtx = matrix_file, 
                       # features = features, 
                      #  cells = barcodes)


#expression_matrix <- ReadSTARsolo(mtx = matrix_file, 
                       # features = features, 
                      #  cells = barcodes)


#seurat_object <- CreateSeuratObject(counts = expression_matrix, project = "DOL_scRNAseq") # Create a Seurat object from the expression matrix and assign it to the variable seurat_object , project name is DOL_scRNAseq

expression_matrix <- ReadSTARsolo(data.dir = "4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered")

seurat_object_filtered <- (counts = expression_matrix, project = "DOL_scRNAseq_filtered") # Create a Seurat object from the expression matrix and assign it to the variable seurat_object , project name is DOL_scRNAseq_filtered

# Sauvegarder l'objet Seurat filtré
saveRDS(seurat_object_filtered, file = "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds")

# Vérifier que le fichier a été créé et afficher des informations
fichier_seurat <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_filtered.rds"

if (file.exists(fichier_seurat)) {
    # Obtenir les informations sur le fichier
    info_fichier <- file.info(fichier_seurat)
    
    # Afficher les informations
    cat("✓ Fichier créé avec succès :", fichier_seurat, "\n")
    cat("Taille du fichier :", round(info_fichier$size/1024/1024, 2), "MB\n")
    cat("Date de création :", format(info_fichier$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Afficher les dimensions de l'objet Seurat
    cat("\nInformations sur l'objet Seurat :\n")
    cat("Nombre de cellules :", ncol(seurat_object_filtered), "\n")
    cat("Nombre de gènes :", nrow(seurat_object_filtered), "\n")
} else {
    cat("❌ Erreur : Le fichier n'a pas été créé.")
}

expression_matrix2 <- ReadSTARsolo(data.dir = "4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/raw")

seurat_object_raw <- CreateSeuratObject(counts = expression_matrix2, project = "DOL_scRNAseq_raw") # Create a Seurat object from the expression matrix and assign it to the variable seurat_object , project name is DOL_scRNAseq_raw

# Sauvegarder l'objet Seurat
saveRDS(seurat_object_raw, file = "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds")

fichier_seurat2 <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq_raw.rds"

if (file.exists(fichier_seurat2)) {
    # Obtenir les informations sur le fichier
    info_fichier2 <- file.info(fichier_seurat2)
    
    # Afficher les informations
    cat("✓ Fichier créé avec succès :", fichier_seurat2, "\n")
    cat("Taille du fichier :", round(info_fichier2$size/1024/1024, 2), "MB\n")
    cat("Date de création :", format(info_fichier2$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Afficher les dimensions de l'objet Seurat
    cat("\nInformations sur l'objet Seurat :\n")
    cat("Nombre de cellules :", ncol(seurat_object_raw), "\n")
    cat("Nombre de gènes :", nrow(seurat_object_raw), "\n")
} else {
    cat("❌ Erreur : Le fichier n'a pas été créé.")
}




####voir avec Solène pour les erreurs de lecture des fichiers mtx UniqueAndMult-Uniform.mtx qui contiennent des float et non des int

