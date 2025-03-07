# Chargement des librairies

    library(reticulate) # Permet d'utiliser des fonctions python dans R
    library(Seurat) # Permet de manipuler des données de single-cell RNA-seq
    library(ggplot2) # Permet de faire des graphiques
    library(dplyr) # Permet de manipuler des données
    library(Matrix) # Permet de manipuler des matrices
    library(cowplot) # Permet de faire des graphiques
    library(forcats) # Permet de manipuler des facteurs
    library(RColorBrewer) # Permet de manipuler les couleurs




#script pour convertir les noms de gènes en noms de gènes utiles pour les bactéries
source("7_Article/script/script_Val/usefull_convert_gene_names_bacteria.R")
#prend en entrée un tableau csv de corespondance entre les noms de gènes et les noms de gènes utiles pour les bactéries ainsi qu'un vecteur de noms de gènes qui est déja fourni par les autres codes. 

#liste les fonctions importées de : 7_Article/script/script_Val/usefull_convert_gene_names_bacteria.R


#script pour assigner les barcodes (cellules) aux puits de la plaque 
source("7_Article/script/script_Val/usefull_assign_cell_wells.R" )
#ce script devra etre modifié si jamais les barcodes ne sont pas dans l'ordre des puits de la plaque dans une autre expérience


#script pour filtrer les barcodes (cellules) les plus fiables avec un threshold à définir 
source("7_Article/script/script_Val/usefull_round_one_bc_collapse.R" ) 
#prend en entrée un csv et data 


#script pour assigner les métadonnées aux cellules pour M14 et M15
source("7_Article/script/script_Val/usefull_assign_metadata.R")

lsf.str()  # listes les fonctions importées 





create_seurat_object_from_seq_files <- function(data_dir, sublibrary, conversion_file_gene_names, csv_file_filter, threshold = NULL, exclude_types, ribosome_file = NULL, min.cells,
    min.features) {
  # TODO: Améliorer la gestion des doublons au lieu de les supprimer
  # TODO: Voir la fonction round_one_bc_collapse pour le seuil de filtrage qui doit etre améliorer pour etre plus simple, et la génération du csv de filtre qui doit etre améliorer pour etre plus simple
  # TODO 

  # FIXME: 

  # Cette fonction crée un objet Seurat à partir de fichiers de données de séquençage

  # data_dir : chemin du répertoire contenant les fichiers de données
  # sublibrary : nom de la sous-bibliothèque
  # threshold : seuil de filtrage des cellules si nécessaire (NULL par défaut)
  # ribosome_file : chemin du fichier de ribosomes
  # conversion_file_gene_names : chemin du fichier de conversion des noms de gènes
  # csv_file_filter : chemin du fichier de filtre des cellules
  

  # Vérifier si le répertoire de données existe
  if (!dir.exists(data_dir)) {
    stop("Le répertoire spécifié n'existe pas : ", data_dir)
  }
  
  # Construire les chemins des fichiers
  matrix_file <- file.path(data_dir, "UniqueAndMult-Uniform.mtx")
  features_file <- file.path(data_dir, "features.tsv")
  barcodes_file <- file.path(data_dir, "barcodes.tsv")

  # Vérifier si les fichiers existent
  if (!file.exists(matrix_file) || !file.exists(features_file) || !file.exists(barcodes_file)) {
    stop("Un ou plusieurs fichiers de données sont introuvables dans : ", data_dir)
  }

  #Verifier si possède les fonctions utilisées : convert_gene_names_bacteria, round_one_bc_collapse, assign_cell_wells
  if (!exists("convert_gene_names_bacteria")) {
    stop("La fonction convert_gene_names_bacteria n'est pas définie.")
  }
  if (!exists("round_one_bc_collapse")) {
    stop("La fonction round_one_bc_collapse n'est pas définie.")
  }
  if (!exists("assign_cell_wells")) {
    stop("La fonction assign_cell_wells n'est pas définie.")
  }
  
  # Vérifier si les packages nécessaires sont installés
  required_packages <- c("Seurat", "Matrix")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Le package ", pkg, " n'est pas installé.")
    }
  })


  # Lecture des fichiers de comptage
  data <- readMM(matrix_file)
  genes <- read.table(features_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
  genes <- convert_gene_names_bacteria(genes, conversion_file_gene_names)
  barcodes <- read.table(barcodes_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
  
  # Structuration des données
  rownames(data) <- genes
  colnames(data) <- barcodes

  # Filtrage des cellules uniquement si un seuil est spécifié
  if (!is.null(threshold)) {
    filteredData <- round_one_bc_collapse(data, threshold, csv_file_filter)
  } else {
    filteredData <- data  # Aucune suppression si threshold est NULL
  }

  # Vérifier s'il y a des doublons et les supprimer
  if (any(duplicated(rownames(filteredData)))) {
    warning("Des gènes en double ont été trouvés et seront supprimés.")
    duplicated_genes <- rownames(filteredData)[duplicated(rownames(filteredData))]
    print(duplicated_genes)
    filteredData <- filteredData[!duplicated(rownames(filteredData)), ]
  }





  # Suppression des types de gènes spécifiés
  if (!is.null(exclude_types) && length(exclude_types) > 0) {
    if (is.null(ribosome_file) || !file.exists(ribosome_file)) {
      stop("Le fichier de classification des gènes est introuvable ou non spécifié : ", ribosome_file)
    }
    gene_classes <- read.csv(ribosome_file, stringsAsFactors = FALSE)
    colnames(gene_classes) <- c('gene', 'type')
    genes_to_exclude <- gene_classes[gene_classes$type %in% exclude_types, "gene"]
    filteredData <- filteredData[!(rownames(filteredData) %in% genes_to_exclude), ]
  }

  # Création de l'objet Seurat
  SO <- CreateSeuratObject(counts = filteredData, 
    #min.cells = 15, #15
    #min.features = 5 #5
  )
# La fonction CreateSeuratObject de Seurat permet de créer un objet Seurat à partir de données brutes. Voici les principaux arguments que vous pouvez utiliser :
# 	•	counts : Un objet de type matrice contenant des données non normalisées, avec les cellules en colonnes et les caractéristiques (gènes) en lignes.
# 	•	assay : Nom de l’assay initial. Par défaut, il est défini sur “RNA”.
# 	•	names.field : Si les noms de vos cellules suivent un format spécifique, cet argument permet de choisir quel champ utiliser pour définir la classe d’identité initiale de chaque cellule. Par exemple, si vos cellules sont nommées comme “BARCODE_CLUSTER_CELLTYPE”, définir names.field à 3 utilisera “CELLTYPE” comme identité initiale.
# 	•	names.delim : Délimiteur utilisé dans les noms de colonnes des cellules pour séparer les différentes parties. Par exemple, si vos cellules sont nommées “BARCODE-CLUSTER-CELLTYPE”, définissez cet argument à “-” pour séparer le nom de la cellule en ses composants.
# 	•	meta.data : Données supplémentaires au niveau des cellules à ajouter à l’objet Seurat. Doit être un data.frame où les lignes correspondent aux noms des cellules et les colonnes aux champs de métadonnées supplémentaires. Les noms de lignes dans les métadonnées doivent correspondre aux noms de colonnes de la matrice counts.
# 	•	project : Nom du projet pour l’objet Seurat.
# 	•	min.cells : Inclut les caractéristiques détectées dans au moins ce nombre de cellules. Cela sous-ensemble également la matrice counts. Pour réintroduire des caractéristiques exclues, créez un nouvel objet avec un seuil inférieur.
# 	•	min.features : Inclut les cellules où au moins ce nombre de caractéristiques est détecté.

# Ces arguments permettent de personnaliser la création de votre objet Seurat en fonction de la structure de vos données et des besoins de votre analyse


# # Ajout de métadonnées
# wells <- assign_cell_wells(colnames(SO@assays$RNA))

# # Assigner directement les données de 'well' et 'sublibrary' au Seurat object
# SO[["well"]] <- wells$well
# SO$sublibrary <- sublibrary

  # Ajout des métadonnées
  wells <- assign_cell_wells(colnames(SO))
  SO <- AddMetaData(SO, metadata = wells, col.name = "well")
  SO <- AddMetaData(SO, metadata = rep(sublibrary, ncol(SO)), col.name = "sublibrary")

  # Affichage des métadonnées
  print("Seurat object créé avec succès")
  print(SO)

  print ("Métadonnées :")

 # Affichage des métadonnées de l'objet Seurat
print("Affichage des 6 premières lignes des métadonnées :")
print(head(SO@meta.data))  # Affiche les 6 premières lignes des métadonnées

print("Dimensions des métadonnées (nombre de cellules et de colonnes) :")
print(dim(SO@meta.data))  # Vérifie le nombre de cellules et de colonnes

print("Noms des colonnes des métadonnées :")
print(colnames(SO@meta.data))  # Liste les colonnes disponibles dans les métadonnées
return(SO)
}