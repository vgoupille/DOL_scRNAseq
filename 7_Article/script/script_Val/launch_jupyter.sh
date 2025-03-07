#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=jupyter.log

#       SBATCH --gres=gpu:1  # Ajoute un GPU si nécessaire

#affiche la memoire et le cpu et le gpu
echo "Memory: $SLURM_MEM_PER_NODE" > jupyter.log
echo "CPU: $SLURM_CPUS_PER_TASK" > jupyter.log
echo "GPU: $CUDA_VISIBLE_DEVICES" > jupyter.log

#permet de voir les noeuds sur lesquels le job a été lancé
scontrol show job $SLURM_JOB_ID | grep NodeList > jupyter.log

#permet de lancer l'image docker pour analyse single cell et jupyter notebook avec un env R , python ....
singularity exec 6_Containers/shortcakelight.sif jupyternotebook.sh
