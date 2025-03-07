

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



# fichier de ribosomes
path_ribosome_file <- '7_Article/script/utile_bact/bacteria.ribosomes.csv'
# fichier de conversion des noms de gènes
path_conversion_file_gene_names <- '7_Article/script/utile_bact/bacteria_gene_conversion.csv'
# fichier de filtre des barcodes (cellules)
path_csv_file_filter <- '7_Article/script/utile_bact/r6ptorderedbcs.csv'


# Définir le chemin du répertoire
plots_dir <- "7_Article/results/plots/"

# Créer le répertoire si nécessaire
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Définir le chemin du fichier PDF de sortie
pdf_output_path <- paste0(plots_dir, "M15_mRNA_all_plots.pdf")


# Ouvrir un fichier PDF pour capturer tous les graphiques
pdf(pdf_output_path, width = 8, height = 6)

cat("📄 Enregistrement des graphiques dans :", pdf_output_path, "\n")

# Appel de la fonction qui génère l'objet Seurat et produit des graphiques

M15_mRNA_nofilter  <- create_seurat_object_from_seq_files(
    data_dir = replica1, 
    sublibrary ='M15',
    threshold = 0.99, 
    ribosome_file = path_ribosome_file,
    conversion_file_gene_names = path_conversion_file_gene_names, csv_file_filter = path_csv_file_filter,
    exclude_types = c('rRNA', 'tRNA'),
    #min.cells = 15, #min.cell permet de filtrer les cellules qui ont un nombre de gènes supérieur à min.cell
    #min.features = 5 #min.features permet de filtrer les gènes qui sont exprimés dans un nombre de cellules supérieur à min.features
    )

# Fermer le fichier PDF pour sauvegarder les graphiques
dev.off()
cat("✅ Les graphiques ont été sauvegardés dans", pdf_output_path, "\n")

# Sauvegarde de l'objet Seurat
saveRDS(M15_mRNA_nofilter, "7_Article/results/Seurat/M15_mRNA_nofilter.rds")



# Chemin du fichier sauvegardé
file_path <- "7_Article/results/Seurat/M15_mRNA_nofilter.rds"

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




#Add metadata to the cells M15 for OD
print ("Assigning metadata to cells M15: ")
M15_mRNA_nofilter_DO <- assign_metadata_OD_M15(M15_mRNA_nofilter)


#Test d'affichage des données
M15_mRNA_nofilter_DO 
head(M15_mRNA_nofilter_DO@meta.data)
dim(M15_mRNA_nofilter_DO@meta.data)  # Vérifie le nombre de cellules et de colonnes
colnames(M15_mRNA_nofilter_DO@meta.data)  # Liste les colonnes disponibles



# Sauvegarde de l'objet Seurat
saveRDS(M15_mRNA_nofilter_DO, "7_Article/results/Seurat/M15_mRNA_nofilter_DO.rds")

# Chemin du fichier sauvegardé
file_path <- "7_Article/results/Seurat/M15_mRNA_nofilter_DO.rds"

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