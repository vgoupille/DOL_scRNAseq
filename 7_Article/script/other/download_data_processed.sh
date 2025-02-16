#!/bin/bash

#SBATCH --job-name=download_data_processed
#SBATCH --ntasks=1                  # Nombre de tâches
#SBATCH --cpus-per-task=8           # Nombre de cœurs par tâche
#SBATCH --mem=20G                   # Mémoire allouée
#SBATCH --time=02:00:00             # Temps d'exécution maximum (hh:mm:ss)
#SBATCH --output=download_data_processed_%j.log  # Fichier de sortie (le %j est remplacé par l'ID du job SLURM)
#SBATCH --error=download_data_processed_%j.err   # Fichier d'erreur

# Vérifier si le dossier cible existe, sinon le créer
if [ ! -d "7_Article/data" ]; then
  echo "Le dossier 7_Article/data n'existe pas. Création du dossier."
  mkdir -p 7_Article/data
fi

# URL du fichier à télécharger
URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE151940&format=file"
OUTPUT="7_Article/data/GSE151940_RAW.tar"

# Vérifier si l'URL est accessible
if curl --head --silent --fail "$URL" > /dev/null; then
  echo "L'URL est accessible, téléchargement en cours..."
  wget -O $OUTPUT $URL

  # Vérifier si le téléchargement a réussi
  if [ -f "$OUTPUT" ]; then
    echo "Téléchargement réussi et fichier sauvegardé sous $OUTPUT"
  else
    echo "Erreur : le téléchargement n'a pas réussi."
    exit 1
  fi
else
  echo "Erreur : l'URL n'est pas accessible."
  exit 1
fi

#pour décompresser le fichier
tar -xvf $OUTPUT -C 7_Article/data