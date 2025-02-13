#srun --mem 20G --pty bash
#source /local/env/envconda.sh 
#conda activate /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/5_environnements/env_seurat


# Notre repertoire de travail/ Faire un WORKINGDIR
WORKINGDIR="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/"
setwd(WORKINGDIR)
getwd()
 
# Load the Seurat package
library(Seurat)
