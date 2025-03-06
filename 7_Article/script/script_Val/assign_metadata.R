
#Fonction pour associé les metadonnnées (ici la densité optique et le numéro des puits) puis les données de séquençage de B. subtilis

assign_metadata_OD_M15 <- function(raw.genefull) {
  print(sum(duplicated(rownames(raw.genefull@meta.data))))  # Vérifie les doublons avant fusion

  # Vérifier si l'objet `raw.genefull` est un SeuratObject
  if (!inherits(raw.genefull, "Seurat")) {
    stop("Erreur : L'objet fourni n'est pas un objet Seurat.")
  }

  # Vérifier la présence de la colonne "well" dans les métadonnées
  if (!"well" %in% colnames(raw.genefull@meta.data)) {
    raw.genefull@meta.data$well <- rownames(raw.genefull@meta.data)  # Utilise les noms des cellules
  }

  # Définition des conditions basées sur les puits (wells)
  conditions <- list(
    'OD0.5' = paste0('A', 1:6),
    'OD1.0' = paste0('A', 7:12),
    'OD1.3' = paste0('B', 1:6),
    'OD1.6' = paste0('B', 7:12),
    'OD2.8' = paste0('C', 1:6),
    'OD3.6' = paste0('C', 7:12),
    'OD5.3' = paste0('D', 1:6),
    'OD6.0' = paste0('D', 7:12)
  )

  # Assigner la condition à chaque cellule
  raw.genefull@meta.data <- raw.genefull@meta.data %>%
    mutate(cond = case_when(
      well %in% conditions$OD0.5 ~ "OD0.5",
      well %in% conditions$OD1.0 ~ "OD1.0",
      well %in% conditions$OD1.3 ~ "OD1.3",
      well %in% conditions$OD1.6 ~ "OD1.6",
      well %in% conditions$OD2.8 ~ "OD2.8",
      well %in% conditions$OD3.6 ~ "OD3.6",
      well %in% conditions$OD5.3 ~ "OD5.3",
      well %in% conditions$OD6.0 ~ "OD6.0",
      TRUE ~ NA_character_
    ))

  print(sum(duplicated(rownames(raw.genefull@meta.data))))  # Vérifie les doublons après fusion
  return(raw.genefull)
}





assign_metadata_OD_M14 <- function(raw.genefull) {
  print(sum(duplicated(rownames(raw.genefull@meta.data))))  # Vérifie les doublons avant fusion

  # Vérification que `raw.genefull` est un objet Seurat
  if (!inherits(raw.genefull, "Seurat")) {
    stop("Erreur : L'objet fourni n'est pas un objet Seurat.")
  }

  # Vérification et ajout de la colonne `well` si absente
  if (!"well" %in% colnames(raw.genefull@meta.data)) {
    raw.genefull@meta.data$well <- rownames(raw.genefull@meta.data)
  }

  # Définition des conditions avec `paste0`
  conditions <- list(
    'OD0.5' = paste0('A', 1:8),
    'OD1.0' = c(paste0('A', 9:12), paste0('B', 1:4)),
    'OD1.7' = paste0('B', 5:12),
    'OD2.0' = paste0('C', 1:8),
    'OD2.8' = c(paste0('C', 9:12), paste0('D', 1:4)),
    'OD3.2' = paste0('D', 5:12)
  )

  # Assigner les conditions directement aux métadonnées
  raw.genefull@meta.data <- raw.genefull@meta.data %>%
    mutate(cond = case_when(
      well %in% conditions$OD0.5 ~ "OD0.5",
      well %in% conditions$OD1.0 ~ "OD1.0",
      well %in% conditions$OD1.7 ~ "OD1.7",
      well %in% conditions$OD2.0 ~ "OD2.0",
      well %in% conditions$OD2.8 ~ "OD2.8",
      well %in% conditions$OD3.2 ~ "OD3.2",
      TRUE ~ NA_character_
    ))

  print(sum(duplicated(rownames(raw.genefull@meta.data))))  # Vérifie les doublons après fusion
  return(raw.genefull)
}

