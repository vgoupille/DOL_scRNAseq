#!/bin/bash

#SBATCH --job-name=STARsolo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=starsolo_%j.out
#SBATCH --error=starsolo_%j.err

# Activation de l'environnement conda
. /local/env/envconda.sh
conda activate 2_environnements/env_STARsolo

# Paramètres
THREADS=16
GENOME_FASTA="1_data/genome_ref/GCA_030064105.1_ASM3006410v1_genomic.fna"
GFF3_FILE="1_data/genome_ref/genomic.gff"
READ1="1_data/Fastq/microSPLIT-600cells_S1_L001_R1_001.fastq.gz"
READ2="1_data/Fastq/microSPLIT-600cells_S1_L001_R2_001.fastq.gz"
GENOME_DIR="genome_index"
OUTPUT_DIR="starsolo_output"

# Chemins des fichiers de barcodes
BARCODE_R1="1_data/barcodes/Barcode1.txt"
BARCODE_R2R3="1_data/barcodes/Barcode2-3.txt"

# Création des répertoires nécessaires
mkdir -p $GENOME_DIR $OUTPUT_DIR

# Vérification des fichiers d'entrée
echo "$(date) - Vérification des fichiers d'entrée..."
for file in "$GENOME_FASTA" "$GFF3_FILE" "$READ1" "$READ2" "$BARCODE_R1" "$BARCODE_R2R3"; do
    if [ ! -f "$file" ]; then
        echo "ERREUR: Fichier manquant: $file"
        exit 1
    fi
done

# 1. Conversion du GFF3 en GTF
echo "$(date) - Conversion GFF3 vers GTF..."
GTF_FILE="${GFF3_FILE%.gff}.gtf"
gffread $GFF3_FILE -T -v -o $GTF_FILE

# 2. Correction du fichier GTF
echo "$(date) - Correction du fichier GTF..."
GTF_FIXED="${GTF_FILE%.gtf}_fixed.gtf"
awk 'BEGIN{OFS=FS="\t"} {if ($3 == "CDS" || $3 == "transcript") $3 = "exon"; print}' $GTF_FILE > $GTF_FIXED

# 3. Génération de l'index du génome
echo "$(date) - Génération de l'index du génome..."
STAR --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $GENOME_DIR \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $GTF_FIXED

# 4. Alignement et comptage avec STARsolo
echo "$(date) - Exécution de l'alignement STARsolo..."
STAR --genomeDir $GENOME_DIR \
    --runThreadN $THREADS \
    --readFilesIn $READ1 $READ2 \
    --readFilesCommand gunzip -c \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
    --soloUMIposition 0_0_0_9 \
    --soloCBwhitelist $BARCODE_R2R3 $BARCODE_R2R3 $BARCODE_R1  \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_All \
    --soloFeatures Gene GeneFull \
    --soloMultiMappers Uniform \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 50 \
    --alignSJDBoverhangMin 1000 \
    --alignSJoverhangMin 1000 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$OUTPUT_DIR/"

# 5. Vérification des fichiers de sortie
echo "$(date) - Vérification des fichiers de sortie..."
if [ -f "$OUTPUT_DIR/Log.final.out" ] && [ -f "$OUTPUT_DIR/Solo.out/GeneFull/Raw/UniqueAndMult-Uniform.mtx" ]; then
    echo "Analyse terminée avec succès!"
    echo "Fichiers importants pour l'analyse en aval :"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/barcodes.tsv"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/features.tsv"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/Raw/UniqueAndMult-Uniform.mtx"
else
    echo "ERREUR: Certains fichiers de sortie sont manquants!"
    exit 1
fi

echo "$(date) - Fin du traitement"