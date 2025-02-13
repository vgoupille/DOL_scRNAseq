#! /bin/bash
#SBATCH --job-name=create_env_seurat

#SBATCH --output=starsolo_%j.out
#SBATCH --error=starsolo_%j.err

. /local/env/envconda.sh

conda create -p /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/5_environnements/env_seurat r-seurat==3.0.2 -y


#https://biocontainers.pro/tools/r-seurat

#docker pull quay.io/biocontainers/r-seurat:3.0.2--r36h0357c0b_0
#singularity run https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_1