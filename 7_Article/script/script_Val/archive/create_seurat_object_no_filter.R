

# Chargement des librairies

    library(reticulate) # Permet d'utiliser des fonctions python dans R
    library(Seurat) # Permet de manipuler des données de single-cell RNA-seq
    library(ggplot2) # Permet de faire des graphiques
    library(dplyr) # Permet de manipuler des données
    library(Matrix) # Permet de manipuler des matrices
    library(cowplot) # Permet de faire des graphiques
    library(forcats) # Permet de manipuler des facteurs
    library(RColorBrewer) # Permet de manipuler les couleurs



getwd()
setwd("/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/")
getwd()


#script pour create_seurat_object_from_seq_files
source ("7_Article/script/script_Val/usefull_create_seurat_object_from_seq_files.R")




# Chargement des fichiers utilisés pour l'analyse : 

## Chargement des données STARsolo pour les deux réplicats avec chacun  :

# barcodes.tsv => les codes-barres des cellules
# features.tsv => les gènes

# UniqueAndMult-Uniform.mtx => les comptages de gènes par cellule (ici en format Matrix Market) et contient les allignements uniques et multiples des reads (un read peut être aligné à plusieurs endroits dans le génome)
# matrix.mtx => les comptages de gènes par cellule (ici en format Matrix Market) et contient les allignements uniques des reads (un read est aligné à un seul endroit dans le génome)

#Gene or genefull (ici GeneFull : Résultats similaires à Gene/, mais en incluant **toutes les régions des gènes** et pas seulement les exons)) 
#raw or filtered data (ici raw : pas de préfiltrage StarSolo)

# Replica n°1 : M15 (GEO661)
replica1 <- '7_Article/data/data_osf/GEO661'
# Replica n°2 : M14 (GEO660)
replica2 <- '7_Article/data/data_osf/GEO660'


# fichier de ribosomes
path_ribosome_file <- '7_Article/script/utile_bact/bacteria.ribosomes.csv'
# fichier de conversion des noms de gènes
path_conversion_file_gene_names <- '7_Article/script/utile_bact/bacteria_gene_conversion.csv'
# fichier de filtre des barcodes (cellules)
path_csv_file_filter <- '7_Article/script/utile_bact/r6ptorderedbcs.csv'



M15_mRNA_no_filter  <- create_seurat_object_from_seq_files(
    data_dir = replica1, 
    sublibrary ='M15',
    threshold = 0, 
    ribosome_file = path_ribosome_file,
    conversion_file_gene_names = path_conversion_file_gene_names, csv_file_filter = path_csv_file_filter,
    exclude_types = c('rRNA', 'tRNA'),
    #min.cells = 15,
    #min.features = 5
    )


# assigne les metadonnées de densité optique aux cellules

print ("Assigning metadata to cells M15: ")
M15_mRNA_no_filter <- assign_metadata_bacteria_M15(M15_mRNA_no_filter)


#Test d'affichage des données

M15_mRNA_no_filter 
head(M15_mRNA_no_filter@meta.data)
dim(M15_mRNA_no_filter@meta.data)  # Vérifie le nombre de cellules et de colonnes
colnames(M15_mRNA_no_filter@meta.data)  # Liste les colonnes disponibles



# Sauvegarde des objets Seurat pour les données brutes avec tous les types d'ARN 
saveRDS(M15_mRNA_no_filter, "7_Article/results/Seurat/M15_mRNA_no_filter.rds")




M15_mRNA_no_filter <- "7_Article/results/Seurat/M15_mRNA_no_filter.rds"

if (file.exists(M15_mRNA_no_filter)) {
    # Obtenir les informations sur le fichier
    info_fichier <- file.info(M15_mRNA_no_filter)
    
    # Afficher les informations
    cat("✓ Fichier créé avec succès :", M15_mRNA_no_filter, "\n")
    cat("Taille du fichier :", round(info_fichier$size/1024/1024, 2), "MB\n")
    cat("Date de création :", format(info_fichier$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Afficher les dimensions de l'objet Seurat
    cat("\nInformations sur l'objet Seurat :\n")
    cat("Nombre de cellules :", ncol(seurat_object_raw), "\n")
    cat("Nombre de gènes :", nrow(seurat_object_raw), "\n")
} else {
    cat("❌ Erreur : Le fichier n'a pas été créé.")





M14_mRNA_no_filter  <- create_seurat_object_from_seq_files(
    data_dir = replica2, 
    sublibrary ='M14',
    threshold = 0, 
    ribosome_file = path_ribosome_file,
    conversion_file_gene_names = path_conversion_file_gene_names, csv_file_filter = path_csv_file_filter,
    exclude_types = c('rRNA', 'tRNA'),
    #min.cells = 15,
    #min.features = 5
    )



#repeat for replicate 2
print ("Assigning metadata to cells M14: ")
M14_mRNA_no_filter <- assign_metadata_bacteria_M14(M14_mRNA_no_filter)




M14_mRNA_no_filter 
head(M14_mRNA_no_filter@meta.data)
dim(M14_mRNA_no_filter@meta.data)  # Vérifie le nombre de cellules et de colonnes
colnames(M14_mRNA_no_filter@meta.data)  # Liste les colonnes disponibles




saveRDS(M14_mRNA_no_filter, "7_Article/results/Seurat/M14_mRNA_no_filter.rds")


M14_mRNA_no_filter <- "7_Article/results/Seurat/M14_mRNA_no_filter.rds"

if (file.exists(M14_mRNA_no_filter)) {
    # Obtenir les informations sur le fichier
    info_fichier <- file.info(M14_mRNA_no_filter)
    
    # Afficher les informations
    cat("✓ Fichier créé avec succès :", M14_mRNA_no_filter, "\n")
    cat("Taille du fichier :", round(info_fichier$size/1024/1024, 2), "MB\n")
    cat("Date de création :", format(info_fichier$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Afficher les dimensions de l'objet Seurat
    cat("\nInformations sur l'objet Seurat :\n")
    cat("Nombre de cellules :", ncol(seurat_object_raw), "\n")
    cat("Nombre de gènes :", nrow(seurat_object_raw), "\n")
} else {
    cat("❌ Erreur : Le fichier n'a pas été créé.")