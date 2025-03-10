# Charger les packages nécessaires
library(Seurat)
library(scCustomize)
library (anndata)


#pas utile car anndata fait deja appel a reticulate pour creer un env python et faire tourner scanpy 
#library(reticulate)
# Crée un environnement conda (ou active un environnement existant)
#conda_create("myenv")
#conda_install("myenv", "anndata")
#conda_install("myenv", "scanpy")


#py_config()  # Vérifie que le bon environnement conda est chargé
#use_condaenv("myenv")  # Active l'environnement conda 'myenv'

#scanpy <- import("scanpy")


# Exemple : Charger un objet Seurat
seurat_object <-readRDS("/Users/valentingoupille/Desktop/M14_mRNA_nofilter_DO_edit.rds")  # Ici, pbmc3k est un objet Seurat d'exemple



# Vérifier l'assay actif
DefaultAssay(seurat_object) <- "RNA"

str (seurat_object)

# Vérifier les couches disponibles dans l'assay RNA
seurat_object[["RNA"]]@layers

names(seurat_object[["RNA"]]@layers)

# Accéder aux données de la couche 'counts' de l'assay 'RNA'
LayerData(seurat_object, assay = "RNA", layer = "counts")

# Convertir l'objet Seurat en AnnData
as.anndata(
  x = seurat_object,
  file_path = "/Users/valentingoupille/Desktop/",
  file_name = "M14_mRNA_nofilter_DO_edit.h5ad",
  main_layer = "counts", other_layers = NULL
)



#Lire l'objet AnnData grace au package scanpy qui est importé avec anndata
adata <- read_h5ad("/Users/valentingoupille/Desktop/M14_mRNA_nofilter_DO_edit.h5ad")

#Vérifier que l'objet a bien été lu
adata

# affiche l'objet AnnData
head(adata$X)

head(adata$obs)
head(adata$var)


head(seurat_object@meta.data$orig.ident)
#probleme orig.ident n'est pas dans l'objet anndata donc besoin de corriger cela 


