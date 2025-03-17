#!/bin/bash
#SBATCH --job-name=R_job_singularity           # Nom du job
#SBATCH --output=R_job_%j.out                  # Fichier de sortie
#SBATCH --error=R_job_%j.err                   # Fichier d'erreur
#SBATCH --time=04:00:00                        # Temps max d'exécution
#SBATCH --cpus-per-task=16                     # Nombre de cœurs
#SBATCH --mem=32G                              # Mémoire allouée

# Chargement automatique de singularity
. /local/env/envsingularity-3.8.5.sh


singularity exec 6_Containers/shortcakelight.sif Rscript 7_Article/script/script_Val/3_combine_rep1_and_rep2.R
