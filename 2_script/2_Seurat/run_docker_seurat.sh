#!/bin/bash
#SBATCH --job-name=analyse_seurat_docker
#SBATCH --ntasks=1                  # Nombre de tâches
#SBATCH --cpus-per-task=4           # Nombre de cœurs par tâche
#SBATCH --mem=20G                   # Mémoire allouée
#SBATCH --time=01:00:00             # Temps d'exécution maximum (hh:mm:ss)
#SBATCH --output=analyse_seurat_%j.log  # Fichier de sortie (le %j est remplacé par l'ID du job SLURM)
#SBATCH --error=analyse_seurat_%j.err   # Fichier d'erreur

#pour lancer le container dans le terminal
#singularity shell /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/6_Containers/seurat.sif 
#Rscript 2_script/2_Seurat/seurat.R

singularity exec /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/6_Containers/seurat.sif Rscript 2_script/2_Seurat/seurat.R