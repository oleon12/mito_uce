#!/bin/bash
#SBATCH --job-name=mafft_cds
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/mafft_cds_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/mafft_cds_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/odl7/sturnira_mito

source ~/miniforge3/etc/profile.d/conda.sh
conda activate mt_pipeline

INPUT="results/combined_masked_cds_intersection.fasta"
OUTDIR="results/alignments"
OUTPUT="$OUTDIR/combined_masked_cds_intersection.mafft.fasta"

mkdir -p "$OUTDIR"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input FASTA not found: $INPUT"
    exit 1
fi

mafft --thread ${SLURM_CPUS_PER_TASK:-8} --auto "$INPUT" > "$OUTPUT"

echo "CDS alignment written to: $OUTPUT"
echo "N sequences:" $(grep -c "^>" "$OUTPUT")
