# Notre repertoire de travail/ Faire un WORKINGDIR
WORKINGDIR="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/"
setwd(WORKINGDIR)
getwd()
 

#Important files which come from the STARsolo output are the ‘barcodes.tsv’, ‘features.tsv’ and ‘UniqueAndMult-Uniform.mtx’ files. These files are used to create the Seurat object. The ‘UniqueAndMult-Uniform.mtx’ file contains the gene expression data in matrix format. The ‘barcodes.tsv’ file contains the cell barcodes and the ‘features.tsv’ file contains the gene names.: 

#Refer to three additional output files for further downstream analysis: the barcode.tsv, features.tsv and UniqueAndMult-Uniform.mtx. The first two files are generated in the ‘Solo.out/GeneFull/’ folder. The UniqueAndMult-Uniform.mtx is located in the ‘Solo.out/ GeneFull/Raw/’ folder.

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/filtered/barcodes.tsv

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/filtered/features.tsv

 #/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/Gene/raw/UniqueAndMult-Uniform.mtx


#▲ CRITICAL STEP Do not use the ‘matrix.mtx’ folder inside of the ‘Solo.out/GeneFull/’ folder. This file contains only the data from the uniquely aligned genes, whereas we use the data from both unique and multiple alignments.

# Charger les bibliothèques nécessaires
library(Seurat)

# Définir les chemins des fichiers
barcodes <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/barcodes.tsv"
features <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/features.tsv"
matrix_file <- "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/raw/UniqueAndMult-Uniform.mtx"

# Lire les données
data_matrix <- ReadMATRIX(mtx = matrix_file, 
                          features = features, 
                          barcodes = barcodes)

# Créer l'objet Seurat
seurat_obj <- CreateSeuratObject(counts = data_matrix, project = "DOL_scRNAseq")

# Sauvegarder l'objet Seurat
saveRDS(seurat_obj, file = "/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/4_results/1_script_results/2_resuts_seurat/DOL_scRNAseq.rds")