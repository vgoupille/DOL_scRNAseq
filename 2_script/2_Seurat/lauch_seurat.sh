#!/bin/bash
#SBATCH --job-name=analyse_seurat
#SBATCH --ntasks=1                  # Nombre de tâches
#SBATCH --cpus-per-task=4           # Nombre de cœurs par tâche
#SBATCH --mem=20G                   # Mémoire allouée
#SBATCH --time=01:00:00             # Temps d'exécution maximum (hh:mm:ss)
#SBATCH --output=analyse_seurat_%j.log  # Fichier de sortie (le %j est remplacé par l'ID du job SLURM)
#SBATCH --error=analyse_seurat_%j.err   # Fichier d'erreur


# Activation de l'environnement conda
. /local/env/envconda.sh

conda activate /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/5_environnements/env_seurat

# Exécuter votre script R
Rscript /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/analyse_seurat.R