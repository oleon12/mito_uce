#!/bin/bash
#SBATCH --job-name=extract_masked_cds
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/extract_masked_cds_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/extract_masked_cds_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -uo pipefail

cd /scratch/odl7/sturnira_mito

############################################
# CONFIG
############################################
CONDA_ENV="mt_pipeline"
SAMPLE_LIST="CONFS/sample_list.txt"
BED="results/cds_coords.bed"
OUTDIR="results"

############################################
# ACTIVATE CONDA
############################################
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

############################################
# ABSOLUTE PATHS
############################################
SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
BED=$(realpath "$BED")
OUTDIR=$(realpath -m "$OUTDIR")

FAILED_LIST="$OUTDIR/extract_masked_cds_failed_samples.txt"
: > "$FAILED_LIST"

############################################
# CHECK INPUTS
############################################
if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: sample list not found: $SAMPLE_LIST"
    exit 1
fi

if [[ ! -f "$BED" ]]; then
    echo "ERROR: BED file not found: $BED"
    exit 1
fi

############################################
# LOOP OVER SAMPLES
############################################
while read -r SAMPLE READ1 READ2; do
    [[ -z "$SAMPLE" ]] && continue

    CONS="$OUTDIR/$SAMPLE/consensus_masked/${SAMPLE}_mt_masked.fasta"
    GENEDIR="$OUTDIR/$SAMPLE/genes_masked"
    OUTFA="$GENEDIR/${SAMPLE}_cds_masked.fasta"
    TMPFA="$GENEDIR/${SAMPLE}_cds_masked.fasta.tmp"

    mkdir -p "$GENEDIR"
    rm -f "$TMPFA"

    if [[ -f "$OUTFA" ]]; then
        echo "Masked CDS already exists for $SAMPLE. Skipping."
        continue
    fi

    echo "========================================"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$CONS" ]]; then
        echo "Missing masked consensus: $SAMPLE" | tee -a "$FAILED_LIST"
        continue
    fi

    bedtools getfasta \
        -fi "$CONS" \
        -bed "$BED" \
        -s \
        -fo "$TMPFA" || {
        echo "bedtools getfasta failed: $SAMPLE" | tee -a "$FAILED_LIST"
        rm -f "$TMPFA"
        continue
    }

    # simple sanity check
    nseq=$(grep -c "^>" "$TMPFA" 2>/dev/null || true)

    if [[ "$nseq" -eq 0 ]]; then
        echo "No CDS extracted: $SAMPLE" | tee -a "$FAILED_LIST"
        rm -f "$TMPFA"
        continue
    fi

    mv "$TMPFA" "$OUTFA"
    echo "Finished: $SAMPLE with $nseq CDS entries"

done < "$SAMPLE_LIST"

echo "Done. Failed samples listed in: $FAILED_LIST"
