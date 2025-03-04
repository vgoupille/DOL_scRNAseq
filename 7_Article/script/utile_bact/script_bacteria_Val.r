library(reticulate)
library(dplyr)
library(cowplot)
library(Seurat) ##Note this code was written for a previous version of Seurat (4.1.1)
library(ggplot2)
library(Matrix)
library(forcats)
library(akmedoids)
library(RColorBrewer)





# import the R script with other functions

source("script_stage_preprint/utile_bact/convert_gene_names.R")


source("script_stage_preprint/utile_bact/assign_cell_wells.R")


source("script_stage_preprint/utile_bact/round_one_bc_collapse.R")

















#--bacteria data processing-----------------------------------------------------------------------------

#function to create Seurat Object for downstream analysis from counts, features, barcodes, and other metadata
create_seurat_object_from_seq_files <- function(data_dir, sublibrary, ribosome_removal, threshold){
  
  setwd(list.dirs(data_dir))
  
  data <- readMM('UniqueAndMult-Uniform.mtx')
  genes <- read.table('features.tsv')
  genes <- genes$V1
  genes <- convert_gene_names_bacteria(genes)
  barcodes <- read.table('barcodes.tsv')
  barcodes <- barcodes$V1 # est le nom automatique de la première colonne
  
  rownames(data) <- genes
  colnames(data) <- barcodes
  
  #use filtered cells

  filteredData <- round_one_bc_collapse(data, threshold)
  
  #deal with ribosomal reads
  
  if (ribosome_removal == 'rRNA') {
    ribosomes <- read.csv('C:/Users/lbrettne/Desktop/bacteria.ribosomes.csv')
    colnames(ribosomes) <- c('gene','type')
    rRNA <- ribosomes[which(ribosomes$type == 'rRNA'),]
    mRNA <- filteredData[which(!(rownames(filteredData) %in% rRNA$gene)),]
    
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = mRNA)
  } else {
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = filteredData, min.cells = 15, min.features = 5)
  }
  
  #add metadata
  wells <- assign_cell_wells(colnames(SO@assays$RNA))
  celldata <- wells$well
  celldata <- as.data.frame(celldata)
  row.names(celldata) <- wells$barcode
  
  SO$well <- celldata$celldata
  SO$sublibrary <- sublibrary
  
  return(SO)
}

# 1. Fonction create_seurat_object_from_seq_files()

# Cette fonction prend en entrée :
# 	•	data_dir : Répertoire contenant les fichiers de données de séquençage.
# 	•	sublibrary : Nom d’une sous-bibliothèque (sublibrary), utilisé pour l’annotation des cellules.
# 	•	ribosome_removal : Option pour filtrer les ARN ribosomiques ('none' = ne pas filtrer, 'rRNA' = enlever les gènes ribosomiques).
# 	•	threshold : Seuil utilisé dans round_one_bc_collapse() pour filtrer les cellules (probablement en fonction du nombre de lectures).








### Création des objets Seurat pour les données de séquençage de B. subtilis

#B. subtilis replicate 1 
GEO661 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO661 Solo.out/bacillus_only/GeneFull/raw/','M15','none',0.85) # avec ribosome 
GEO661.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO661 Solo.out/bacillus_only/GeneFull/raw/','M15','rRNA',0.85) # sans ribosome

#B. subtilis replicate 2
GEO660 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO660 Solo.out/Bacillus Only/GeneFull/raw/','M14','none',0.875)
GEO660.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO660 Solo.out/Bacillus Only/GeneFull/raw/','M14','rRNA',0.875)




#Fonction pour associé les metadonnnées (ici la densité optique et le numéro des puits) puis les données de séquençage de B. subtilis




# Les deux fonctions, assign_metadata_bacteria et assign_metadata_bacteria_M14, ont un objectif similaire : ajouter des métadonnées (notamment la densité optique cond) aux cellules de l’objet Seurat en fonction du puits (well) où elles ont été détectées. Cependant, il y a des différences dans la structure des regroupements et les conditions de croissance (cond).

# Différences principales :
# 	1.	Nombre de groupes et leurs compositions :
# 	•	assign_metadata_bacteria divise les puits en 8 groupes (s1 à s8), correspondant aux valeurs OD0.5, OD1.0, OD1.3, OD1.6, OD2.8, OD3.6, OD5.3, et OD6.0.
# 	•	assign_metadata_bacteria_M14 divise les puits en 6 groupes (s1 à s6), avec des valeurs légèrement différentes de OD : OD0.5, OD1.0, OD1.7, OD2.0, OD2.8, et OD3.2.
# 	2.	Répartition des puits :
# 	•	assign_metadata_bacteria attribue chaque condition (cond) à des puits spécifiques en ligne A à D (ex : s1 correspond aux puits A1-A6 pour OD0.5).
# 	•	assign_metadata_bacteria_M14 fait de même mais en incluant plus de puits dans certains groupes (ex : s1 inclut A1-A8 au lieu de A1-A6).
# 	3.	Erreurs potentielles dans assign_metadata_bacteria_M14 :
# 	•	Dans s3, le puits B10 apparaît deux fois dans merge(), ce qui est probablement une erreur.
# 	•	La progression des valeurs OD ne suit pas exactement la même logique que assign_metadata_bacteria.

# Conclusion :
# 	•	assign_metadata_bacteria et assign_metadata_bacteria_M14 effectuent la même tâche mais avec une classification différente des puits et des conditions de croissance (OD).
# 	•	assign_metadata_bacteria_M14 a une organisation différente des puits et contient une erreur (B10 en double).
# 	•	Le choix entre les deux dépend du protocole expérimental et des conditions spécifiques de croissance de Bacillus subtilis.



#functions adds cell metadata based on which well barcodes are detected
assign_metadata_bacteria <- function(raw.genefull) {

  s1 <- merge(subset(raw.genefull, subset = well == 'A1'), y = c(subset(raw.genefull, subset = well == 'A2'),
                                                                 subset(raw.genefull, subset = well == 'A3'),
                                                                 subset(raw.genefull, subset = well == 'A4'),
                                                                 subset(raw.genefull, subset = well == 'A5'),
                                                                 subset(raw.genefull, subset = well == 'A6')))
  s1$cond <- 'OD0.5'
  
  s2 <- merge(subset(raw.genefull, subset = well == 'A7'), y = c(subset(raw.genefull, subset = well == 'A8'),
                                                                 subset(raw.genefull, subset = well == 'A9'),
                                                                 subset(raw.genefull, subset = well == 'A10'),
                                                                 subset(raw.genefull, subset = well == 'A11'),
                                                                 subset(raw.genefull, subset = well == 'A12')))
  s2$cond <- 'OD1.0'
  
  s3 <- merge(subset(raw.genefull, subset = well == 'B1'), y = c(subset(raw.genefull, subset = well == 'B2'),
                                                                 subset(raw.genefull, subset = well == 'B3'),
                                                                 subset(raw.genefull, subset = well == 'B4'),
                                                                 subset(raw.genefull, subset = well == 'B5'),
                                                                 subset(raw.genefull, subset = well == 'B6')))
  s3$cond <- 'OD1.3'
  
  s4 <- merge(subset(raw.genefull, subset = well == 'B7'), y = c(subset(raw.genefull, subset = well == 'B8'),
                                                                 subset(raw.genefull, subset = well == 'B9'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B11'),
                                                                 subset(raw.genefull, subset = well == 'B12')))
  s4$cond <- 'OD1.6'
  
  s5 <- merge(subset(raw.genefull, subset = well == 'C1'), y = c(subset(raw.genefull, subset = well == 'C2'),
                                                                 subset(raw.genefull, subset = well == 'C3'),
                                                                 subset(raw.genefull, subset = well == 'C4'),
                                                                 subset(raw.genefull, subset = well == 'C5'),
                                                                 subset(raw.genefull, subset = well == 'C6')))
  s5$cond <- 'OD2.8'
  
  s6 <- merge(subset(raw.genefull, subset = well == 'C7'), y = c(subset(raw.genefull, subset = well == 'C8'),
                                                                 subset(raw.genefull, subset = well == 'C9'),
                                                                 subset(raw.genefull, subset = well == 'C10'),
                                                                 subset(raw.genefull, subset = well == 'C11'),
                                                                 subset(raw.genefull, subset = well == 'C12')))
  s6$cond <- 'OD3.6'
  
  s7 <- merge(subset(raw.genefull, subset = well == 'D1'), y = c(subset(raw.genefull, subset = well == 'D2'),
                                                                 subset(raw.genefull, subset = well == 'D3'),
                                                                 subset(raw.genefull, subset = well == 'D4'),
                                                                 subset(raw.genefull, subset = well == 'D5'),
                                                                 subset(raw.genefull, subset = well == 'D6')))
  s7$cond <- 'OD5.3'
  
  s8 <- merge(subset(raw.genefull, subset = well == 'D7'), y = c(subset(raw.genefull, subset = well == 'D8'),
                                                                 subset(raw.genefull, subset = well == 'D9'),
                                                                 subset(raw.genefull, subset = well == 'D10'),
                                                                 subset(raw.genefull, subset = well == 'D11'),
                                                                 subset(raw.genefull, subset = well == 'D12')))
  s8$cond <- 'OD6.0'
  
  growcurv <- merge(s1,y = c(s2,s3,s4,s5,s6,s7,s8))
  
  return(growcurv)
}



assign_metadata_bacteria_M14 <- function(raw.genefull) {
  
  s1 <- merge(subset(raw.genefull, subset = well == 'A1'), y = c(subset(raw.genefull, subset = well == 'A2'),
                                                                 subset(raw.genefull, subset = well == 'A3'),
                                                                 subset(raw.genefull, subset = well == 'A4'),
                                                                 subset(raw.genefull, subset = well == 'A5'),
                                                                 subset(raw.genefull, subset = well == 'A6'),
                                                                 subset(raw.genefull, subset = well == 'A7'),
                                                                 subset(raw.genefull, subset = well == 'A8')))
  s1$cond <- 'OD0.5'
  
  s2 <- merge(subset(raw.genefull, subset = well == 'A9'), y = c(subset(raw.genefull, subset = well == 'A10'),
                                                                 subset(raw.genefull, subset = well == 'A11'),
                                                                 subset(raw.genefull, subset = well == 'A12'),
                                                                 subset(raw.genefull, subset = well == 'B1'),
                                                                 subset(raw.genefull, subset = well == 'B2'),
                                                                 subset(raw.genefull, subset = well == 'B3'),
                                                                 subset(raw.genefull, subset = well == 'B4')))
  s2$cond <- 'OD1.0'
  
  s3 <- merge(subset(raw.genefull, subset = well == 'B5'), y = c(subset(raw.genefull, subset = well == 'B6'),
                                                                 subset(raw.genefull, subset = well == 'B7'),
                                                                 subset(raw.genefull, subset = well == 'B8'),
                                                                 subset(raw.genefull, subset = well == 'B9'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B12')))
  s3$cond <- 'OD1.7'
  
  s4 <- merge(subset(raw.genefull, subset = well == 'C1'), y = c(subset(raw.genefull, subset = well == 'C2'),
                                                                 subset(raw.genefull, subset = well == 'C3'),
                                                                 subset(raw.genefull, subset = well == 'C4'),
                                                                 subset(raw.genefull, subset = well == 'C5'),
                                                                 subset(raw.genefull, subset = well == 'C6'),
                                                                 subset(raw.genefull, subset = well == 'C7'),
                                                                 subset(raw.genefull, subset = well == 'C8')))
  s4$cond <- 'OD2.0'
  
  s5 <- merge(subset(raw.genefull, subset = well == 'C9'), y = c(subset(raw.genefull, subset = well == 'C10'),
                                                                 subset(raw.genefull, subset = well == 'C11'),
                                                                 subset(raw.genefull, subset = well == 'C12'),
                                                                 subset(raw.genefull, subset = well == 'D1'),
                                                                 subset(raw.genefull, subset = well == 'D2'),
                                                                 subset(raw.genefull, subset = well == 'D3'),
                                                                 subset(raw.genefull, subset = well == 'D4')))
  s5$cond <- 'OD2.8'
  
  s6 <- merge(subset(raw.genefull, subset = well == 'D5'), y = c(subset(raw.genefull, subset = well == 'D6'),
                                                                 subset(raw.genefull, subset = well == 'D7'),
                                                                 subset(raw.genefull, subset = well == 'D8'),
                                                                 subset(raw.genefull, subset = well == 'D9'),
                                                                 subset(raw.genefull, subset = well == 'D10'),
                                                                 subset(raw.genefull, subset = well == 'D11'),
                                                                 subset(raw.genefull, subset = well == 'D12')))
  s6$cond <- 'OD3.2'
  
  growcurv <- merge(s1,y = c(s2,s3,s4,s5,s6))
  return(growcurv)
}

growcurvb <- assign_metadata_bacteria(GEO661)
growcurvb.nr <- assign_metadata_bacteria(GEO661.nr)

#subset data into individual timepoints
OD0.5 <- subset(growcurvb, subset = cond == "OD0.5")
OD1.0 <- subset(growcurvb, subset = cond == "OD1.0")
OD1.3 <- subset(growcurvb, subset = cond == "OD1.3")
OD1.6 <- subset(growcurvb, subset = cond == "OD1.6")
OD2.8 <- subset(growcurvb, subset = cond == "OD2.8")
OD3.6 <- subset(growcurvb, subset = cond == "OD3.6")
OD5.3 <- subset(growcurvb, subset = cond == "OD5.3")
OD6.0 <- subset(growcurvb, subset = cond == "OD6.0")

#subset data minus rRNA into individual timepoints
OD0.5.nr <- subset(growcurvb.nr, subset = cond == "OD0.5")
OD1.0.nr <- subset(growcurvb.nr, subset = cond == "OD1.0")
OD1.3.nr <- subset(growcurvb.nr, subset = cond == "OD1.3")
OD1.6.nr <- subset(growcurvb.nr, subset = cond == "OD1.6")
OD2.8.nr <- subset(growcurvb.nr, subset = cond == "OD2.8")
OD3.6.nr <- subset(growcurvb.nr, subset = cond == "OD3.6")
OD5.3.nr <- subset(growcurvb.nr, subset = cond == "OD5.3")
OD6.0.nr <- subset(growcurvb.nr, subset = cond == "OD6.0")












#fonction pour normaliser 


#normalize and scale mRNA datasets
scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}


# La fonction scale_normalize_seurat_data() effectue les étapes suivantes sur un objet Seurat :
# 	1.	Normalisation des données → NormalizeData()
# 	•	Transforme les valeurs d’expression des gènes en CPM (Counts Per Million).
# 	•	Applique une normalisation logarithmique (LogNormalize par défaut).
# 	2.	Identification des gènes les plus variables → FindVariableFeatures()
# 	•	Sélectionne les nvariableFeatures gènes les plus variables avec la méthode vst (Variance Stabilizing Transformation).
# 	3.	Mise à l’échelle des données → ScaleData()
# 	•	Centre chaque gène sur une moyenne de 0.
# 	•	Met à l’échelle les valeurs d’expression avec un écart-type de 1.
# 	4.	Analyse en Composantes Principales (PCA) → RunPCA()
# 	•	Effectue une PCA sur npcs composantes principales.
# 	•	Utilise les features spécifiés ou, si features = NULL, les gènes variables détectés.
# 	5.	Affichage d’un Elbow Plot → ElbowPlot()
# 	•	Génère un Elbow Plot pour visualiser le nombre optimal de composantes à conserver.
# 	6.	Retourne l’objet Seurat modifié contenant les résultats de la normalisation, du scaling et de la PCA.

# 💡 Résumé :
# 👉 Normalise les données → Sélectionne les gènes variables → Met à l’échelle → Effectue une PCA → Affiche un Elbow Plot.









#Utilisation de la focnrtion scale_normalize_seurat_data pour normaliser et mettre à l'échelle les données de séquençage de B. subtilis

OD0.5.nr.int <- scale_normalize_seurat_data(OD0.5.nr, 3000, 20, NULL)
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.6.nr.int <- scale_normalize_seurat_data(OD1.6.nr, 3000, 20, NULL)
OD2.8.nr.int <- scale_normalize_seurat_data(OD2.8.nr, 3000, 20, NULL)
OD3.6.nr.int <- scale_normalize_seurat_data(OD3.6.nr, 3000, 20, NULL)
OD5.3.nr.int <- scale_normalize_seurat_data(OD5.3.nr, 3000, 20, NULL)
OD6.0.nr.int <- scale_normalize_seurat_data(OD6.0.nr, 3000, 20, NULL)


















#calculate counts and statistics for RNA types, including rRNA, mRNA, etc.
#calculate counts and statistics for RNA types, including rRNA, mRNA, etc.
assign_RNA_types_b <- function(SO, SO2, group_name, timepoint) {
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/bacteria.ribosomes.csv')
  colnames(ribosomes) = c("name","type")
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  tRNA.genes <- ribosomes[which(ribosomes$type == 'tRNA'),]
  rProt.genes <- ribosomes[which(ribosomes$type == 'rProtein'),]

  class <- read.csv('C:/Users/lbrettne/Desktop/bacteria_RP_iESR_RiBi.csv')
  RP.genes <- class[which(class$class == 'RP'),]
  iESR.genes <- class[which(class$class == 'iESR'),]
  RiBi.genes <- class[which(class$class == 'RiBi'),]
  
  RP  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RP.genes$gene),])
  rpcv <- sd(RP[which(RP>0)])/mean(RP[which(RP>0)])
  iESR  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% iESR.genes$gene),])
  RiBi  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RiBi.genes$gene),])
  totalmRNA <- colSums(SO2@assays$RNA)
  
  SO2$rsect <- 0
  SO2$psect <- 0
  
  for (i in 1:dim(SO2)[2]) {
    
    SO2$rsect[i] <- (RP[i]+RiBi[i])/(RP[i]+RiBi[i]+iESR[i])
    SO2$psect[i] <- (iESR[i])/(RP[i]+RiBi[i]+iESR[i])
    
    
  }
  
  dt <- c(0,1.3509309,1.2972221,1.2699011,1.2351942,1.2072780,1.1572952,1.1383160,1.0740787,0.9988778,0.9714859,0.8519719,0.7547571,0.6883084,0.3290760,0.1606523)
  
  rRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rRNA.genes$name),])
  tRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% tRNA.genes$name),])
  mRNA  <- SO@assays$RNA[which(!(rownames(SO@assays$RNA) %in% rRNA.genes$name)),]
  #mRNA  <- SO[which(!(rownames(mRNA) %in% tRNA.genes$name)),]
  mRNA  <- colSums(mRNA)
  rProt <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rProt.genes$name),])
  rProt.ratio <- rProt/mRNA
  totalRNA <- colSums(SO@assays$RNA)
  rRNA.ratio <- rRNA/totalRNA
  mRNA.ratio <- mRNA/totalRNA
  mRNA.rRNA <- mRNA/rRNA
  bulkrr <- sum(rRNA)/sum(totalRNA)

  
  bulkrp <- log10(sum(rProt)/sum(totalRNA)*10^6)
  rcv <- sd(rRNA)/mean(rRNA)
  rrcv <- sd(rRNA.ratio)/mean(rRNA.ratio)
  #rpcv <- sd(rProt[which(rProt>0)])/mean(rProt[which(rProt>0)])
  rr.med.dev.med <- median(abs(rRNA - median(rRNA)))
  
  temp <- data.frame(group_name = group_name, 
                     rRNA = rRNA, 
                     tRNA = tRNA, 
                     mRNA = mRNA, 
                     rProt = rProt,
                     totalRNA = totalRNA,
                     rProt.ratio = rProt.ratio,
                     rRNA.ratio = rRNA.ratio,
                     mRNA.ratio = mRNA.ratio,
                     mRNA.rRNA = mRNA.rRNA,
                     bulkrr = bulkrr,
                     #GpH = 1/(dt[timepoint]),
                     GpH = dt[timepoint],
                     bulkrp = bulkrp,
                     rcv = rcv,
                     rpcv = rpcv,
                     rrcv = rrcv,
                     rsd = sd(rRNA),
                     rrsd = sd(rRNA.ratio),
                     rpsd = sd(rProt),
                     rprsd = sd(rProt.ratio),
                     rrm = mean(rRNA),
                     rr.med.dev.med = rr.med.dev.med,
                     rsect = SO2$rsect,
                     psect = SO2$psect)
  
  return(temp)
}

# La fonction assign_RNA_types_b() effectue une analyse des différents types d’ARN présents dans un objet Seurat (SO2) et calcule des statistiques associées. Voici les étapes détaillées :

# 🔹 1. Lecture des fichiers CSV contenant les classifications des ARN
# 	•	Fichier bacteria.ribosomes.csv : Contient une liste de gènes associés aux types d’ARN ribosomiques (rRNA, tRNA, protéines ribosomiques).
# 	•	Stocke les gènes rRNA, tRNA et protéines ribosomiques (rProtein).
# 	•	Fichier bacteria_RP_iESR_RiBi.csv : Contient des classes fonctionnelles des gènes.
# 	•	Sélectionne les gènes appartenant aux catégories :
# 	•	RP (Ribosomal Proteins, protéines ribosomiques).
# 	•	iESR (Immediate Environmental Stress Response, réponse au stress environnemental).
# 	•	RiBi (Ribosome Biogenesis, biogenèse du ribosome).

# 🔹 2. Calcul des sommes d’expression des différents types d’ARN
# 	•	Pour chaque type d’ARN, on somme les comptes d’expression dans l’objet Seurat SO2@assays$RNA :
# 	•	RP → Protéines ribosomiques
# 	•	iESR → Réponse rapide au stress
# 	•	RiBi → Biogenèse des ribosomes
# 	•	totalmRNA → Somme totale de l’expression de tous les ARNs
# 	•	Calcul du coefficient de variation des protéines ribosomiques (rpcv) :
# 	•	rpcv = sd(RP[which(RP>0)]) / mean(RP[which(RP>0)])
# 	•	Permet d’évaluer la variabilité des protéines ribosomiques parmi les cellules analysées.

# 🔹 3. Calcul des ratios rsect et psect pour chaque cellule

# Ces deux scores évaluent l’équilibre entre les différents types d’ARN dans chaque cellule :
# 	•	rsect (Ratio de l’expression ribosomique sur l’ensemble ARN total) :
# ￼
# → Plus la valeur est élevée, plus la cellule est engagée dans la synthèse ribosomique.
# 	•	psect (Ratio d’ARN liés au stress sur l’ensemble ARN total) :
# ￼
# → Plus la valeur est élevée, plus la cellule est engagée dans une réponse au stress.

# Ces valeurs sont calculées pour chaque cellule de SO2.

# 🔹 Résumé des fonctionnalités

# ✅ Lit des fichiers CSV contenant les gènes associés aux protéines ribosomiques, ARN de stress et biogenèse du ribosome.
# ✅ Calcule l’expression totale des différents types d’ARN.
# ✅ Évalue la variabilité de l’expression des protéines ribosomiques (rpcv).
# ✅ Assigne à chaque cellule des scores (rsect, psect) indiquant si elle est plus orientée vers la production ribosomique ou la réponse au stress.

# 🧐 Interprétation
# 	•	Cellules avec un rsect élevé → Fortement engagées dans la production de ribosomes.
# 	•	Cellules avec un psect élevé → En réponse au stress environnemental.

# 💡 Pourquoi c’est utile ?
# 👉 Cette fonction est précieuse pour analyser l’état transcriptionnel des cellules en fonction de leur engagement dans la croissance active (synthèse ribosomique) ou la réponse au stress.





# La fonction assign_RNA_types_b calcule différentes statistiques et totaux pour des types d’ARN spécifiques (comme les rRNA, mRNA, tRNA, etc.) à partir de données d’un objet Seurat (SO) et d’un second objet Seurat (SO2). Voici un résumé détaillé de ce qu’elle fait :

# 1. Lecture des fichiers CSV :
# 	•	La fonction commence par lire deux fichiers CSV qui contiennent des informations sur les types de gènes (ribosomes, protéines ribosomales) et leurs classes (RP, iESR, RiBi), respectivement dans bacteria.ribosomes.csv et bacteria_RP_iESR_RiBi.csv.
# 	•	Ensuite, elle sépare les gènes en fonction de leur type ou classe (rRNA, tRNA, rProtein, RP, iESR, RiBi).

# 2. Calcul des totaux pour les gènes spécifiques :
# 	•	Pour chaque type de gène (RP, iESR, RiBi), la fonction utilise colSums() pour additionner les valeurs d’expression des gènes correspondant dans l’objet SO2.
# 	•	Elle calcule également le total de tous les gènes (totalmRNA).

# 3. Calcul des ratios de sections :
# 	•	Les variables rsect et psect sont calculées pour chaque échantillon à l’aide de ces formules :
# 	•	rsect[i] = (RP + RiBi) / (RP + RiBi + iESR)
# 	•	psect[i] = iESR / (RP + RiBi + iESR)
# 	•	Ces ratios reflètent la proportion des différentes classes de gènes dans chaque échantillon.

# 4. Calcul des autres types d’ARN :
# 	•	Ensuite, elle calcule les totaux pour les gènes associés aux rRNA, tRNA et mRNA en utilisant colSums() pour chaque catégorie de gènes (rRNA, tRNA, rProtein).
# 	•	Les ratios suivants sont calculés pour chaque échantillon :
# 	•	rProt.ratio = rProtéine / mRNA
# 	•	rRNA.ratio = rRNA / totalRNA
# 	•	mRNA.ratio = mRNA / totalRNA
# 	•	mRNA.rRNA = mRNA / rRNA

# 5. Calculs supplémentaires :
# 	•	bulkrr : La proportion globale de rRNA dans l’ARN total.
# 	•	bulkrp : Le log10 de la proportion de protéines ribosomales dans l’ARN total.
# 	•	Diverses autres statistiques sont calculées, comme les coefficients de variation (rcv), les médianes et écarts-types pour les gènes spécifiques.

# 6. Création de la dataframe finale :
# 	•	La fonction crée une dataframe avec toutes les statistiques calculées pour chaque échantillon, en ajoutant aussi des informations sur le groupe (group_name) et les différentes métriques (ratios, écarts-types, médianes, etc.).
# 	•	Elle retourne cette dataframe avec toutes les informations calculées pour un groupe donné à un point temporel spécifique.

# 7. Paramètres :
# 	•	SO: L’objet Seurat principal contenant les données d’expression des gènes (ARN).
# 	•	SO2: Un autre objet Seurat qui est utilisé pour calculer des statistiques supplémentaires basées sur d’autres types de gènes (par exemple, RP, iESR, RiBi).
# 	•	group_name: Le nom du groupe pour lequel les calculs sont effectués.
# 	•	timepoint: Le point temporel pour lequel le calcul de la croissance des cellules est réalisé.

# Résumé global :

# La fonction calcule des statistiques détaillées pour les types d’ARN (rRNA, tRNA, mRNA, etc.), les totaux d’expressions, les ratios spécifiques (par exemple, rProt/mRNA), et des indices de variabilité. Elle utilise ces données pour créer un tableau des métriques pour chaque échantillon à un temps donné et le renvoie sous forme de dataframe.

# Est-ce que cela répond à ta question sur le fonctionnement de cette fonction ? 😊











  


OD0.5.RNA <- assign_RNA_types_b(OD0.5, OD0.5.nr.int, "0.5", 4)
OD1.0.RNA <- assign_RNA_types_b(OD1.0, OD1.0.nr.int, "1.0", 5)
OD1.3.RNA <- assign_RNA_types_b(OD1.3, OD1.3.nr.int, "1.3", 8)
OD1.6.RNA <- assign_RNA_types_b(OD1.6, OD1.6.nr.int, "1.6", 9)
OD2.8.RNA <- assign_RNA_types_b(OD2.8, OD2.8.nr.int, "2.8", 12)
OD3.6.RNA <- assign_RNA_types_b(OD3.6, OD3.6.nr.int, "3.6", 14)
OD5.3.RNA <- assign_RNA_types_b(OD5.3, OD5.3.nr.int, "5.3", 15)
OD6.0.RNA <- assign_RNA_types_b(OD6.0, OD6.0.nr.int, "6.0", 16)

odRNAb <- rbind(OD0.5.RNA,OD1.0.RNA,OD1.3.RNA,OD1.6.RNA,OD2.8.RNA,OD3.6.RNA,OD5.3.RNA,OD6.0.RNA)







#repeat for replicate 2
growcurvb.m14 <- assign_metadata_bacteria_M14(GEO660)
growcurvb.nr.m14 <- assign_metadata_bacteria_M14(GEO660.nr)

OD0.5.m14 <- subset(growcurvb.m14, subset = cond == "OD0.5")
OD1.0.m14 <- subset(growcurvb.m14, subset = cond == "OD1.0")
OD1.7.m14 <- subset(growcurvb.m14, subset = cond == "OD1.7")
OD2.0.m14 <- subset(growcurvb.m14, subset = cond == "OD2.0")
OD2.8.m14 <- subset(growcurvb.m14, subset = cond == "OD2.8")
OD3.2.m14 <- subset(growcurvb.m14, subset = cond == "OD3.2")

OD0.5.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD0.5")
OD1.0.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD1.0")
OD1.7.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD1.7")
OD2.0.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD2.0")
OD2.8.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD2.8")
OD3.2.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD3.2")

OD0.5.nr.m14.int <- scale_normalize_seurat_data(OD0.5.nr.m14, 3000, 20, NULL)
OD1.0.nr.m14.int <- scale_normalize_seurat_data(OD1.0.nr.m14, 3000, 20, NULL)
OD1.7.nr.m14.int <- scale_normalize_seurat_data(OD1.7.nr.m14, 3000, 20, NULL)
OD2.0.nr.m14.int <- scale_normalize_seurat_data(OD2.0.nr.m14, 3000, 20, NULL)
OD2.8.nr.m14.int <- scale_normalize_seurat_data(OD2.8.nr.m14, 3000, 20, NULL)
OD3.2.nr.m14.int <- scale_normalize_seurat_data(OD3.2.nr.m14, 3000, 20, NULL)

m14OD0.5.RNA <- assign_RNA_types_b(OD0.5.m14, OD0.5.nr.m14.int, "0.5", 4)
m14OD1.0.RNA <- assign_RNA_types_b(OD1.0.m14, OD1.0.nr.m14.int, "1.0", 5)
m14OD1.7.RNA <- assign_RNA_types_b(OD1.7.m14, OD1.7.nr.m14.int, "1.7", 10)
m14OD2.0.RNA <- assign_RNA_types_b(OD2.0.m14, OD2.0.nr.m14.int, "2.0", 11)
m14OD2.8.RNA <- assign_RNA_types_b(OD2.8.m14, OD2.8.nr.m14.int, "2.8", 12)
m14OD3.2.RNA <- assign_RNA_types_b(OD3.2.m14, OD3.2.nr.m14.int, "3.2", 13)

odRNAb.m14 <- rbind(m14OD0.5.RNA,m14OD1.0.RNA,m14OD1.7.RNA,m14OD2.0.RNA,m14OD2.8.RNA,m14OD3.2.RNA)





#####

#--Figure 1-----------------------------------------------------------------------------



doublingrate <- function(fit, x, t){
  drs <- c()
  for (i in 2:length(fit)){
    dr <- log2(fit[i]/fit[i-1])/(x[i]-x[i-1])
    drs <- c(drs,dr)
  }
  
  grs <- c()
  for (j in 2:length(t)){
    k <- which.min(abs(x - t[j]))
    gr <- drs[k]
    grs <- c(grs, gr)
  }
  return(grs)
}



#bacteria
hour  <- c(0.00,1.10,2.13,2.45,2.75,2.90,3.18,3.28,3.60,3.93,4.0,4.42,4.67,4.95,6.12,6.97) #Kuchina and Brettner et al. 2021 (Science)
OD600 <- c(0.02,0.04,0.27,0.49,0.94,0.93,1.16,1.34,1.58,1.76,2.0,2.87,3.20,3.56,5.31,6.34)

fitb <- nls(OD600 ~ L/(1+exp(-k*(hour - x0))), start = list(L = 6, k = 0.5, x0 = 3))
Lb  <- summary(fitb)$parameters[1,1]
kb  <- summary(fitb)$parameters[2,1]
x0b <- summary(fitb)$parameters[3,1]
xb <- seq(from = 0, to = 8, by = 0.1)
fitb <- Lb/(1+exp(-kb*(xb - x0b)))

GpHb <- doublingrate(fitb,xb,xb)
GpHb <- data.frame(t = xb[-1], gph = GpHb)
GpHb_data <- doublingrate(fitb,xb,hour)
GpHb_data <- data.frame(t = hour[c(4,5,8,9,12,13,14,15)], gph = GpHb_data[c(3,4,7,8,11,12,13,14)])

cell_density <- data.frame(t = xb, cpml = fitb)
cell_dens_data <- data.frame(t = hour[c(4,5,8,9,12,13,14,15)], cpml = OD600[c(4,5,8,9,12,13,14,15)] )

b1.1.1 <- ggplot(GpHb, aes(x = t, y = gph)) +
  geom_line() +
  geom_point(data = GpHb_data, aes(x = t, y = gph), size = 4, fill = "#8CA2CB", shape = 21) +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  scale_y_continuous(position = "right") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.1.1

b1.1.2 <- ggplot(cell_density, aes(x = t, y = cpml)) +
  geom_line(color = "gray") +
  geom_point(data = cell_dens_data, aes(x = t, y = cpml), size = 4, fill = "#8CA2CB", shape = 21, color = "gray") +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.1.2

b1.2 <- ggplot(odRNAb, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.2

x <- c(mean(OD0.5.RNA$GpH), 
       mean(OD1.0.RNA$GpH), 
       mean(OD1.3.RNA$GpH), 
       mean(OD1.6.RNA$GpH),
       mean(OD2.8.RNA$GpH),
       mean(OD3.6.RNA$GpH),
       mean(OD5.3.RNA$GpH),
       mean(OD6.0.RNA$GpH))

y <- c(mean(OD0.5.RNA$rRNA), 
       mean(OD1.0.RNA$rRNA), 
       mean(OD1.3.RNA$rRNA), 
       mean(OD1.6.RNA$rRNA),
       mean(OD2.8.RNA$rRNA),
       mean(OD3.6.RNA$rRNA),
       mean(OD5.3.RNA$rRNA),
       mean(OD6.0.RNA$rRNA))

summary(lm(y~x))




b1.3 <- ggplot(odRNAb, aes(x = GpH, y = log10(rRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(2,3,4,5),labels = c("100","1,000", "10,000", "100,000")) +
  NoLegend()
b1.3

summary(lm(log10(odRNAb$rRNA)~odRNAb$GpH))



# ce code génère plusieurs graphiques à l’aide de ggplot2. Voici un résumé des graphiques créés dans le code :

# 	1.	Graphique b1.1.1 :
# 	•	Titre : Taux de croissance en fonction du temps.
# 	•	Ce graphique montre le taux de croissance (générations par heure, gph) calculé par la fonction doublingrate en fonction du temps (t), avec une courbe lissée de ce taux de croissance. Des points sont ajoutés pour les valeurs observées à certains moments spécifiques (heure hour).

# 	2.	Graphique b1.1.2 :
# 	•	Titre : Densité cellulaire en fonction du temps.
# 	•	Ce graphique montre la densité cellulaire (mesurée par OD600), qui représente l’absorbance à 600 nm, en fonction du temps (t). Les valeurs ajustées à partir du modèle sigmoïde sont tracées, avec des points indiquant les valeurs réelles mesurées à certains moments du temps.

# 	3.	Graphique b1.2 :
# 	•	Titre : Relation entre le taux de croissance et la quantité d’ARN ribosomique par cellule.
# 	•	Ce graphique examine la relation entre le taux de croissance (en générant GpH) et la quantité d’ARN ribosomique (rRNA) par cellule. Les points moyens (en fonction de GpH) sont visualisés et un modèle linéaire est ajusté pour montrer cette relation.

# 	4.	Graphique b1.3 :
# 	•	Titre : Distribution de log10(rRNA) en fonction du taux de croissance.
# 	•	Ce graphique montre la distribution de la quantité d’ARN ribosomique en utilisant log10(rRNA) en fonction du taux de croissance (GpH). Il combine un graphique de type violon (distribution des données) avec des points de résumé pour la moyenne et une ligne de régression linéaire ajustée.

# Ces graphiques sont générés de manière transparente avec un fond et des axes stylisés, et certains ont des ajustements spécifiques comme l’absence de légende ou une échelle log pour l’axe des ordonnées.




















#--Figure 2-----------------------------------------------------------------------------

###======> doit conrespondre à la figure 3 c de l'article 

#1. Graphique b2.1 : Variation du coefficient de variation des ARN ribosomiques


b2.1 <- ggplot(odRNAb, aes(x = GpH, y = rcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  #stat_summary(data = odRNAb.m14, aes(x = GpH, y = rcv, fill = factor(GpH)), fun = mean, geom = "point", size = 4, fill = "#2c4470", shape = 21) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend() +
  xlab("population GR") +
  ylab("rRNA counts/cell CV")
b2.1




b2.2 <- ggplot(odRNAb, aes(x = GpH, y = rpcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  #stat_summary(data = odRNAb.m14, aes(x = GpH, y = rcv, fill = factor(GpH)), fun = mean, geom = "point", size = 4, fill = "#2c4470", shape = 21) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
  ) +
  NoLegend() +
  ylab("rProtein mRNA CV") +
  xlab("population GR")
b2.2


# La différence entre les deux graphiques (b2.1 et b2.2) réside principalement dans les variables mesurées sur l’axe des y, bien que la structure générale du graphique soit identique. Voici les détails de chaque graphique :

# Graphique b2.1 :
# 	•	Axe des y : Le graphique montre la variation du coefficient de variation (CV) des comptages d’ARN ribosomiques par cellule (rcv). Le rcv est une mesure qui reflète la variation relative des comptages d’ARN ribosomiques dans les cellules.
# 	•	Axe des x : Il s’agit du taux de croissance de la population en générations par heure (GpH).
# 	•	But : Ce graphique examine si la variation des comptages d’ARN ribosomiques (rRNA) au niveau cellulaire est liée au taux de croissance de la population. Le coefficient de variation des rRNA par cellule est représenté en fonction du taux de croissance de la population.

# Graphique b2.2 :
# 	•	Axe des y : Ce graphique montre la variation du coefficient de variation (CV) des mRNA des protéines ribosomiques par cellule (rpcv). Le rpcv est une mesure similaire, mais elle concerne spécifiquement les mRNA des protéines ribosomiques (par opposition aux ARN ribosomiques dans le graphique précédent).
# 	•	Axe des x : Comme dans b2.1, l’axe des x représente le taux de croissance de la population en générations par heure (GpH).
# 	•	But : Ce graphique examine si la variation des mRNA des protéines ribosomiques au niveau cellulaire est liée au taux de croissance de la population.

# Résumé des différences :
# 	1.	Variable mesurée en y :
# 	•	b2.1 mesure le coefficient de variation des comptages d’ARN ribosomiques (rRNA).
# 	•	b2.2 mesure le coefficient de variation des mRNA des protéines ribosomiques (rProtein mRNA).
# 	2.	Interprétation :
# 	•	Dans b2.1, vous explorez la variation des rRNA dans les cellules en fonction du taux de croissance de la population.
# 	•	Dans b2.2, vous explorez la variation des mRNA des protéines ribosomiques dans les cellules, en fonction également du taux de croissance de la population.

# Structure similaire :
# 	•	Les deux graphiques ont la même structure avec des points moyens (calculés via stat_summary(fun = mean)), une ligne de régression linéaire (geom_smooth(method = "lm")) et un fond transparent.





#--Figure 3---------------------------------------------------------------------



scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}
find_neighbors_clusters <- function(integratedSO, dims, resolution){
  
  integrated <- FindNeighbors(integratedSO, dims = 1:dims)
  clustered <- FindClusters(integrated, resolution = resolution)
  clustered <- RunUMAP(clustered, dims = 1:dims)
  
  return(clustered)
}
differentially_expressed_genes <- function(clusteredSO) {
  
  markers <- FindAllMarkers(clusteredSO, only.pos = TRUE)
  
  for (j in 0:(length(levels(markers$cluster))-1))
  {
    cluster.genes <- markers[which(markers$cluster == j),]
    cluster.genes <- cluster.genes[which(cluster.genes$p_val_adj < 0.05),]
    print(head(cluster.genes[order(-cluster.genes$avg_log2FC),], n = 100))
  }
}

###ici fait graph de cluster UMAP 
#bacteria 
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.0.nr.clust <- find_neighbors_clusters(OD1.0.nr.int, 10, 0.1)
OD1.0.nr.clust$rRNA <- OD1.0.RNA$rRNA

OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.3.nr.clust <- find_neighbors_clusters(OD1.3.nr.int, 10, 0.1)
OD1.3.nr.clust$rRNA <- OD1.3.RNA$rRNA

b3.1 <- DimPlot(OD1.0.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#339966", "#808080"))
b3.1

b3.2 <- DimPlot(OD1.3.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#339966", "#808080"))
b3.2

OD1.0.zero <- subset(OD1.0.nr.clust, subset = seurat_clusters == 0)
OD1.0.one <- subset(OD1.0.nr.clust, subset = seurat_clusters == 1)

OD1.0.01 <- rbind(data.frame(rRNA = OD1.0.zero$rRNA, cluster = 0), 
                 data.frame(rRNA = OD1.0.one$rRNA, cluster = 1))

OD1.3.zero <- subset(OD1.3.nr.clust, subset = seurat_clusters == 0)
OD1.3.one <- subset(OD1.3.nr.clust, subset = seurat_clusters == 1)

OD1.3.01 <- rbind(data.frame(rRNA = OD1.3.zero$rRNA, cluster = 0), 
                  data.frame(rRNA = OD1.3.one$rRNA, cluster = 1))

b3.2 <- ggplot(OD1.0.01, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
  geom_jitter(alpha = 0.05) +
  #geom_boxplot(notch = TRUE, alpha = 0.8) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_fill_manual(values = c("#339966", "#808080")) +
  theme_classic() +
  ylab("rRNA counts/cell") +
  xlab("cluster") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
b3.2

b3.4 <- ggplot(OD1.3.01, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
  geom_jitter(alpha = 0.05) +
  #geom_boxplot(notch = TRUE, alpha = 0.8) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_fill_manual(values = c("#339966", "#808080")) +
  theme_classic() +
  ylab("rRNA counts/cell") +
  xlab("cluster") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
b3.4

t.test(log10(OD1.0.zero$rRNA),log10(OD1.0.one$rRNA))

differentially_expressed_genes(OD1.0.nr.clust)
differentially_expressed_genes(OD1.3.nr.clust)





#--Figure 4---------------------------------------------------------------------

rRNAquarts <- function(odRNA){
  GpH <- levels(as.factor(odRNA$GpH))
  
  mins <- c()
  maxs <- c()
  
  for (i in 1:length(GpH)) {
    rRNA <- odRNA$rRNA[which(odRNA$GpH == GpH[i])]
    
    mins[i] <- min(rRNA)
    maxs[i] <- max(rRNA)
  }
  
  maxmins = max(mins)
  minmaxs = min(maxs)
  
  rRNArange <- odRNA$rRNA[which(odRNA$rRNA >= maxmins & odRNA$rRNA <= minmaxs)]
  names.r1 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[1] & odRNA$rRNA < summary(rRNArange)[2]),]), quartile = 1)
  names.r2 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[2] & odRNA$rRNA < summary(rRNArange)[3]),]), quartile = 2)
  names.r3 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[3] & odRNA$rRNA < summary(rRNArange)[5]),]), quartile = 3)
  names.r4 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[5] & odRNA$rRNA < summary(rRNArange)[6]),]), quartile = 4)
  
  quartile_bcs <- rbind(names.r1, names.r2, names.r3, names.r4)
  
  return(quartile_bcs)
}







#bacteria




quartiles <- rRNAquarts(odRNAb)

bQ1 <- growcurvb.nr[,which(colnames(growcurvb.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 1)])]
bQ4 <- growcurvb.nr[,which(colnames(growcurvb.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 4)])]

bQ1.int <- scale_normalize_seurat_data(bQ1, 3000, 20, NULL)
bQ4.int <- scale_normalize_seurat_data(bQ4, 3000, 20, NULL)

bQ1.clustered <- find_neighbors_clusters(bQ1.int, 10, 0.2)
bQ4.clustered <- find_neighbors_clusters(bQ4.int, 10, 0.1)

b4.1 <- DimPlot(bQ4.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#99CCB2", "#808080"))
b4.1

#differentially_expressed_genes(bQ4.clustered)

b4.2 <-DimPlot(bQ4.clustered, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666","#90E0EF","#B47EDE")) + 
  labs(title = NULL)
b4.2

b4.3 <- DimPlot(bQ1.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#808080","#99CCB2"))
b4.3

#differentially_expressed_genes(bQ1.clustered)

b4.4 <-DimPlot(bQ1.clustered, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666","#90E0EF","#B47EDE")) + 
  labs(title = NULL)
b4.4

plot_grid(y4.1,b4.1,y4.3,b4.3)
plot_grid(y4.2,b4.2,y4.4,b4.4)














#--Figure 5---------------------------------------------------------------------

b5.1 <- ggplot(test, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b5.1

b5.2 <- ggplot(test, aes(x = GpH, y = log10(mRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(2,3,4,5),labels = c("100","1,000", "10,000", "100,000")) +
  NoLegend()
b5.2

growcurvb.nr.m14.int <- scale_normalize_seurat_data(growcurvb.nr.m14, 3000, 20, NULL)
growcurvb.nr.m14.clust <- find_neighbors_clusters(growcurvb.nr.m14.int, 10, 0.05)
growcurvb.nr.m14.clust$mRNA.rRNA <- odRNAb.m14$mRNA.rRNA

b5.3 <- DimPlot(growcurvb.nr.m14.clust, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 0.01) + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666"))+ NoLegend()
b5.3

DimPlot(growcurvb.nr.m14.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1)

b5.4 <- FeaturePlot(growcurvb.nr.m14.clust, reduction = "umap", features = "mRNA.rRNA", cols = c("black","#8CA2CB", "#89BB87"), pt.size = 0.01)+ 
  theme(aspect.ratio = 1) + NoLegend()
b5.4

b5.3+b5.4

differentially_expressed_genes(growcurvb.nr.m14.clust)

spores <- subset(growcurvb.nr.m14.clust, subset = seurat_clusters == 2)
nsspores <- subset(growcurvb.nr.m14.clust, subset = seurat_clusters != 2)

test0 <- subset(nine.nr.clust, subset = seurat_clusters == 0) 
test1 <- subset(nine.nr.clust, subset = seurat_clusters == 1)

test <- merge(test0, y = test1)
test.int <- scale_normalize_seurat_data(test, 3000, 20, NULL)
test.clust <- find_neighbors_clusters(test.int, 10, 0.2)

differentially_expressed_genes(test.clust)





#--Supplemental Figures---------------------------------------------------------
#--SFs additional means and distributions-------------------------------------------------------------------------




#bacteria

growcurvb.int <- scale_normalize_seurat_data(growcurvb, 3000, 20, NULL)

OD0.5.int <- scale_normalize_seurat_data(OD0.5, 3000, 20, NULL)
OD1.0.int <- scale_normalize_seurat_data(OD1.0, 3000, 20, NULL)
OD1.3.int <- scale_normalize_seurat_data(OD1.3, 3000, 20, NULL)
OD1.6.int <- scale_normalize_seurat_data(OD1.6, 3000, 20, NULL)
OD2.8.int <- scale_normalize_seurat_data(OD2.8, 3000, 20, NULL)
OD3.6.int <- scale_normalize_seurat_data(OD3.6, 3000, 20, NULL)
OD5.3.int <- scale_normalize_seurat_data(OD5.3, 3000, 20, NULL)
OD6.0.int <- scale_normalize_seurat_data(OD6.0, 3000, 20, NULL)

OD0.5.RNA.int <- assign_RNA_types_b(OD0.5.int, "0.5", 4)
OD1.0.RNA.int <- assign_RNA_types_b(OD1.0.int, "1.0", 5)
OD1.3.RNA.int <- assign_RNA_types_b(OD1.3.int, "1.3", 8)
OD1.6.RNA.int <- assign_RNA_types_b(OD1.6.int, "1.6", 9)
OD2.8.RNA.int <- assign_RNA_types_b(OD2.8.int, "2.8", 12)
OD3.6.RNA.int <- assign_RNA_types_b(OD3.6.int, "3.6", 13)
OD5.3.RNA.int <- assign_RNA_types_b(OD5.3.int, "5.3", 14)
OD6.0.RNA.int <- assign_RNA_types_b(OD6.0.int, "6.0", 15)

odRNAb.int <- rbind(OD0.5.RNA.int,OD1.0.RNA.int,OD1.3.RNA.int,OD1.6.RNA.int,OD2.8.RNA.int,OD3.6.RNA.int,OD5.3.RNA.int,OD6.0.RNA.int)

bs1.1 <- ggplot(odRNAb.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.1

x <- c(mean(OD0.5.RNA.int$GpH), 
       mean(OD1.0.RNA.int$GpH), 
       mean(OD1.3.RNA.int$GpH), 
       mean(OD1.6.RNA.int$GpH),
       mean(OD2.8.RNA.int$GpH),
       mean(OD3.6.RNA.int$GpH),
       mean(OD5.3.RNA.int$GpH),
       mean(OD6.0.RNA.int$GpH))

y <- c(mean(OD0.5.RNA.int$rRNA), 
       mean(OD1.0.RNA.int$rRNA), 
       mean(OD1.3.RNA.int$rRNA), 
       mean(OD1.6.RNA.int$rRNA),
       mean(OD2.8.RNA.int$rRNA),
       mean(OD3.6.RNA.int$rRNA),
       mean(OD5.3.RNA.int$rRNA),
       mean(OD6.0.RNA.int$rRNA))

summary(lm(y~x))

bs1.2 <- ggplot(odRNAb.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.2

summary(lm(odRNAb.int$rRNA~odRNAb.int$GpH))

bs1.3 <- ggplot(odRNAb, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.3

x <- c(mean(OD0.5.RNA$GpH), 
       mean(OD1.0.RNA$GpH), 
       mean(OD1.3.RNA$GpH), 
       mean(OD1.6.RNA$GpH),
       mean(OD2.8.RNA$GpH),
       mean(OD3.6.RNA$GpH),
       mean(OD5.3.RNA$GpH),
       mean(OD6.0.RNA$GpH))

y <- c(mean(OD0.5.RNA$mRNA.rRNA), 
       mean(OD1.0.RNA$mRNA.rRNA), 
       mean(OD1.3.RNA$mRNA.rRNA), 
       mean(OD1.6.RNA$mRNA.rRNA),
       mean(OD2.8.RNA$mRNA.rRNA),
       mean(OD3.6.RNA$mRNA.rRNA),
       mean(OD5.3.RNA$mRNA.rRNA),
       mean(OD6.0.RNA$mRNA.rRNA))

summary(lm(y~x))

bs1.4 <- ggplot(odRNAb, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.4

summary(lm(odRNAb$mRNA.rRNA~odRNAb$GpH))

OD0.5.nr.int <- scale_normalize_seurat_data(OD0.5.nr, 3000, 20, NULL)
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.6.nr.int <- scale_normalize_seurat_data(OD1.6.nr, 3000, 20, NULL)
OD2.8.nr.int <- scale_normalize_seurat_data(OD2.8.nr, 3000, 20, NULL)
OD3.6.nr.int <- scale_normalize_seurat_data(OD3.6.nr, 3000, 20, NULL)
OD5.3.nr.int <- scale_normalize_seurat_data(OD5.3.nr, 3000, 20, NULL)
OD6.0.nr.int <- scale_normalize_seurat_data(OD6.0.nr, 3000, 20, NULL)

OD0.5.RNA.nr.int <- assign_RNA_types_b(OD0.5.nr.int, "0.5", 4)
OD1.0.RNA.nr.int <- assign_RNA_types_b(OD1.0.nr.int, "1.0", 5)
OD1.3.RNA.nr.int <- assign_RNA_types_b(OD1.3.nr.int, "1.3", 8)
OD1.6.RNA.nr.int <- assign_RNA_types_b(OD1.6.nr.int, "1.6", 9)
OD2.8.RNA.nr.int <- assign_RNA_types_b(OD2.8.nr.int, "2.8", 12)
OD3.6.RNA.nr.int <- assign_RNA_types_b(OD3.6.nr.int, "3.6", 13)
OD5.3.RNA.nr.int <- assign_RNA_types_b(OD5.3.nr.int, "5.3", 14)
OD6.0.RNA.nr.int <- assign_RNA_types_b(OD6.0.nr.int, "6.0", 15)

odRNAb.nr.int <- rbind(OD0.5.RNA.nr.int,
                       OD1.0.RNA.nr.int,
                       OD1.3.RNA.nr.int,
                       OD1.6.RNA.nr.int,
                       OD2.8.RNA.nr.int,
                       OD3.6.RNA.nr.int,
                       OD5.3.RNA.nr.int,
                       OD6.0.RNA.nr.int)

bs1.5 <- ggplot(odRNAb.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.5

x <- c(mean(OD0.5.RNA.nr.int$GpH), 
       mean(OD1.0.RNA.nr.int$GpH), 
       mean(OD1.3.RNA.nr.int$GpH), 
       mean(OD1.6.RNA.nr.int$GpH),
       mean(OD2.8.RNA.nr.int$GpH),
       mean(OD3.6.RNA.nr.int$GpH),
       mean(OD5.3.RNA.nr.int$GpH),
       mean(OD6.0.RNA.nr.int$GpH))

y <- c(mean(OD0.5.RNA.nr.int$rProt), 
       mean(OD1.0.RNA.nr.int$rProt), 
       mean(OD1.3.RNA.nr.int$rProt), 
       mean(OD1.6.RNA.nr.int$rProt),
       mean(OD2.8.RNA.nr.int$rProt),
       mean(OD3.6.RNA.nr.int$rProt),
       mean(OD5.3.RNA.nr.int$rProt),
       mean(OD6.0.RNA.nr.int$rProt))

summary(lm(y~x))

bs1.6 <- ggplot(odRNAb.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.6

summary(lm(odRNAb.nr.int$rProt~odRNAb.nr.int$GpH))