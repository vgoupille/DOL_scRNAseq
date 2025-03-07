#!/bin/bash
#SBATCH --job-name=R_job_singularity           # Nom du job
#SBATCH --output=R_job_%j.out       # Fichier de sortie (avec ID du job)
#SBATCH --error=R_job_%j.err        # Fichier d'erreur (avec ID du job)
#SBATCH --time=04:00:00             # Temps max d'exécution
#SBATCH --cpus-per-task=16           # Nombre de cœurs
#SBATCH --mem=32G                    # Mémoire allouée


#Chargement automatique de singularity
. /local/env/envsingularity-3.8.5.sh

# Exécuter le script R dans le conteneur 

#script.R est le script R à exécuter de test 
singularity exec 6_Containers/shortcakelight.sif Rscript 7_Article/script/script_Val/script.R


