#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --ntasks=1                  # Nombre de tâches
#SBATCH --cpus-per-task=8           # Nombre de cœurs par tâche
#SBATCH --mem=20G                   # Mémoire allouée
#SBATCH --time=02:00:00             # Temps d'exécution maximum (hh:mm:ss)
#SBATCH --output=download_sra_%j.log  # Fichier de sortie (le %j est remplacé par l'ID du job SLURM)
#SBATCH --error=download_sra_%j.err   # Fichier d'erreur


# rendre le script exécutable
# chmod +x 7_Article/download.sh


# Charger l'environnement SRA Toolkit (ajuste si nécessaire)
. /local/env/envsra-tools-2.11.0.sh

# Définir le répertoire de destination
DEST_DIR="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/7_Article/data"

# Vérifier si le dossier existe, sinon le créer
if [ ! -d "$DEST_DIR" ]; then
    echo "Le dossier n'existe pas. Création de $DEST_DIR..."
    mkdir -p "$DEST_DIR"
fi

# Se déplacer dans le dossier de destination
cd "$DEST_DIR" || exit

# Liste des identifiants SRA
SRA_IDS=("SRR11940660" )
#"SRR11940661" "SRR11940662")



# Télécharger et convertir les fichiers SRA
echo "Téléchargement et conversion des fichiers SRA dans $DEST_DIR..."
for SRA in "${SRA_IDS[@]}"; do
    echo "Téléchargement de $SRA..."
    prefetch "$SRA"

    echo "Conversion de $SRA en FASTQ..."
    fastq-dump --split-files --gzip --outdir "$DEST_DIR" "$SRA"

    # Vérifier si la conversion a réussi
    if [ $? -eq 0 ]; then
        echo "Suppression du fichier .sra de $SRA pour économiser de l'espace."
        rm "$SRA.sra"
    else
        echo "⚠ Erreur lors de la conversion de $SRA ! Vérifiez les logs."
    fi
done

echo "Téléchargement et conversion terminés dans $DEST_DIR."