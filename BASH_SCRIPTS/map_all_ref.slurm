#!/bin/bash
#SBATCH --job-name=map_all_ref
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=244G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/map_all_ref_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/map_all_ref_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

# ---- safety ----
set -euo pipefail

# ---- move to scratch ----
cd /scratch/odl7/sturnira_mito

############################################
# CONFIG
############################################
CONDA_ENV="mt_pipeline"
SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
OUTDIR="results"
THREADS=${SLURM_CPUS_PER_TASK:-8}

############################################
# ACTIVATE CONDA
############################################
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

############################################
# ABSOLUTE PATHS
############################################
SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
OUTDIR=$(realpath -m "$OUTDIR")

############################################
# CHECK INPUTS
############################################
if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: sample list not found: $SAMPLE_LIST"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference FASTA not found: $REF"
    exit 1
fi

mkdir -p "$OUTDIR"

############################################
# LOOP OVER SAMPLES
############################################
while read -r SAMPLE READ1 READ2; do

    [[ -z "$SAMPLE" ]] && continue

    READ1=$(realpath "$READ1")
    READ2=$(realpath "$READ2")

    if [[ ! -f "$READ1" ]]; then
        echo "ERROR: R1 not found for $SAMPLE -> $READ1"
        exit 1
    fi

    if [[ ! -f "$READ2" ]]; then
        echo "ERROR: R2 not found for $SAMPLE -> $READ2"
        exit 1
    fi

    echo "========================================"
    echo "Processing sample: $SAMPLE"
    echo "R1: $READ1"
    echo "R2: $READ2"

    mkdir -p "$OUTDIR/$SAMPLE/bam"

    BAM="$OUTDIR/$SAMPLE/bam/$SAMPLE.sorted.bam"
    MAPPED="$OUTDIR/$SAMPLE/bam/$SAMPLE.mapped.bam"

    if [[ -f "$MAPPED" && -f "${MAPPED}.bai" ]]; then
        echo "Mapped BAM already exists for $SAMPLE. Skipping."
        continue
    fi

    bwa mem -t "$THREADS" "$REF" "$READ1" "$READ2" \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$BAM" -

    samtools index "$BAM"

    samtools view -b -F 4 "$BAM" > "$MAPPED"
    samtools index "$MAPPED"

    echo "Finished: $SAMPLE"

done < "$SAMPLE_LIST"

echo "ALL SAMPLES FINISHED"

