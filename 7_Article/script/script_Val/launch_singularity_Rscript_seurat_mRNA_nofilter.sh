#!/bin/bash
#SBATCH --job-name=R_job_singularity           # Nom du job
#SBATCH --output=R_job_%j.out                  # Fichier de sortie
#SBATCH --error=R_job_%j.err                   # Fichier d'erreur
#SBATCH --time=04:00:00                        # Temps max d'exécution
#SBATCH --cpus-per-task=16                     # Nombre de cœurs
#SBATCH --mem=32G                              # Mémoire allouée

# Chargement automatique de singularity
. /local/env/envsingularity-3.8.5.sh

# Soumettre le premier job SLURM pour le premier script
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=R_job_M15
#SBATCH --output=R_job_M15_%j.out
#SBATCH --error=R_job_M15_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
singularity exec 6_Containers/shortcakelight.sif Rscript 7_Article/script/script_Val/create_seurat_object_nofilter_mRNA_M15.R
EOF

# Soumettre le deuxième job SLURM pour le deuxième script
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=R_job_M14
#SBATCH --output=R_job_M14_%j.out
#SBATCH --error=R_job_M14_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
singularity exec 6_Containers/shortcakelight.sif Rscript 7_Article/script/script_Val/create_seurat_object_nofilter_mRNA_M14.R
EOF