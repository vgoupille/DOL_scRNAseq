# Utiliser renv pour garantir les versions des packages
install.packages("renv")
library(renv)

# Récupérer l'environnement figé si déjà existant
if (file.exists("renv.lock")) {
    # Activer la gestion des versions avec renv
  renv::restore()
} else {
  # Installer des packages CRAN avec versions précises
  renv::install("tidyverse@2.0.0")
  renv::install("data.table@1.14.8")
  renv::install("devtools@2.4.5")
  renv::install("BiocManager@1.30.22")

  # Installer des packages Bioconductor avec versions
  BiocManager::install("DESeq2", version = "3.16")
  BiocManager::install("edgeR", version = "3.42")
  BiocManager::install("limma", version = "3.56")

  # Installer des packages depuis GitHub avec versions
  devtools::install_github("satijalab/seurat", ref = "v4.3.0")  # Version spécifique
  devtools::install_github("thomasp85/ggforce")  # Depuis le repo officiel
  
  
  # Sauvegarder l'état de l'environnement
  renv::snapshot()
}

# Installer IRkernel pour Jupyter
install.packages("IRkernel")

# Ajouter le kernel R à Jupyter (pour qu'il apparaisse dans Jupyter Lab)
IRkernel::installspec(user = FALSE, name = "myenv_R", displayname = "R (myenv_R)")