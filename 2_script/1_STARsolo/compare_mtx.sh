#!/bin/bash

# Vérification des arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 fichier1.mtx fichier2.mtx"
    exit 1
fi

file1="$1"
file2="$2"

# Nettoyage des fichiers (suppression des lignes de commentaires et d'en-tête)
clean_file() {
    grep -v '^%' "$1" | tail -n +2 | sort
}

clean_file "$file1" > file1_clean.mtx
clean_file "$file2" > file2_clean.mtx

# Comparaison des fichiers ligne par ligne
echo "📌 Comparaison des fichiers $file1 et $file2"

# Trouver les différences
diff_output=$(diff file1_clean.mtx file2_clean.mtx)

if [ -z "$diff_output" ]; then
    echo "✅ Les fichiers sont identiques."
else
    echo "❌ Différences trouvées :"
    echo "$diff_output"
fi

# Nettoyage des fichiers temporaires
rm file1_clean.mtx file2_clean.mtx

#sbatch 2_script/1_STARsolo/compare_mtx.sh 4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/matrix.mtx 4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/raw/matrix.mtx


#compare le nombre de cellules dans chaque matrice
#wc -l 4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/raw/barcodes.tsv
#wc -l 4_results/1_script_results/1_results_STARsolo/starsolo_output/Solo.out/GeneFull/filtered/barcodes.tsv