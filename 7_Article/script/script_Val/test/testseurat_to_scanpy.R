# Charger les packages nécessaires
library(Seurat)
library(scCustomize)
library (anndata)


#library(reticulate)
# Crée un environnement conda (ou active un environnement existant)
#conda_create("myenv")
#conda_install("myenv", "anndata")


#py_config()  # Vérifie que le bon environnement conda est chargé
#use_condaenv("myenv")  # Active l'environnement conda 'myenv'

# Exemple : Charger un objet Seurat
seurat_object <-readRDS("/Users/valentingoupille/Desktop/Combined_M14_M15_mRNA_nofilter_DO.rds")  # Ici, pbmc3k est un objet Seurat d'exemple

# Convertir l'objet Seurat en AnnData
as.anndata(x = seurat_object, file_path = "~/Desktop", file_name = "/Users/valentingoupille/Desktop/Combined_M14_M15_mRNA.h5ad")

# Convertir l'objet Seurat en AnnData
as.anndata(
  x = seurat_object,
  file_path = "/Users/valentingoupille/Desktop/",
  file_name = "Combined_M14_M15_mRNA.h5ad",
  main_layer = "combined_counts", other_layers = NULL
)


# Liste des assays disponibles dans ton objet Seurat
names(seurat_object)


# Vérification des cellules dans l'assay
head(seurat_object[["combined_counts"]]@cell.names)

# Vérification des caractéristiques dans l'assay
head(seurat_object[["combined_counts"]]@feature.names)



validObject(seurat_object)



# Vérification des couches disponibles dans l'assay
names(seurat_object[["RNA"]]@layers)

# Vérification des cellules dans l'assay RNA
head(seurat_object[["RNA"]]@cells)

# Vérification des caractéristiques dans l'assay RNA
head(seurat_object[["RNA"]]@features)




# Récupérer les noms des cellules
cell_names <- colnames(seurat_object[["RNA"]]@counts)

# Récupérer les noms des caractéristiques
feature_names <- rownames(seurat_object[["RNA"]]@counts)

# Vérifier les premières valeurs
head(cell_names)
head(feature_names)




#V5 


# Accéder aux couches de comptage
counts_1 <- seurat_object[["RNA"]]@layers$counts.1
counts_2 <- seurat_object[["RNA"]]@layers$counts.2
combined_counts <- seurat_object[["RNA"]]@layers$combined_counts

# Récupérer les noms des cellules et des caractéristiques
cell_names <- colnames(combined_counts)
feature_names <- rownames(combined_counts)

# Vérifier les premiers éléments
head(cell_names)
head(feature_names)



LayerData(seurat_object, assay = "RNA", layer = "combined_counts")

