#!/bin/bash
#SBATCH --job-name=prep_ref
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=244G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/prep_ref_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/prep_ref_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail


cd /scratch/odl7/sturnira_mito

############################################
# CONFIG
############################################
CONDA_ENV="mt_pipeline"
REF="references/S_ludovici_QCAZ_18312.fasta"

############################################
# ACTIVATE CONDA
############################################
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

############################################
# ABSOLUTE PATH
############################################
REF=$(realpath "$REF")

############################################
# CHECK INPUT
############################################
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF"
    exit 1
fi

echo "Reference FASTA: $REF"

############################################
# BUILD INDEXES
############################################
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Building BWA index..."
    bwa index "$REF"
else
    echo "BWA index already exists."
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "Building samtools FASTA index..."
    samtools faidx "$REF"
else
    echo "FASTA index already exists."
fi

echo "Reference preparation complete."

