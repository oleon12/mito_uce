#!/bin/bash
#SBATCH --job-name=cons_mt_all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/cons_mt_all_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/cons_mt_all_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail

cd /scratch/odl7/sturnira_mito

CONDA_ENV="mt_pipeline"
SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
OUTDIR="results"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
OUTDIR=$(realpath -m "$OUTDIR")

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: sample list not found: $SAMPLE_LIST"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference FASTA not found: $REF"
    exit 1
fi

while read -r SAMPLE READ1 READ2; do
    [[ -z "$SAMPLE" ]] && continue

    VCF="$OUTDIR/$SAMPLE/vcf/$SAMPLE.vcf.gz"
    CONSDIR="$OUTDIR/$SAMPLE/consensus"
    CONS="$CONSDIR/${SAMPLE}_mt.fasta"

    if [[ ! -f "$VCF" ]]; then
        echo "WARNING: VCF not found for $SAMPLE -> $VCF"
        continue
    fi

    mkdir -p "$CONSDIR"

    if [[ -f "$CONS" ]]; then
        echo "Consensus already exists for $SAMPLE. Skipping."
        continue
    fi

    echo "========================================"
    echo "Processing sample: $SAMPLE"
    echo "VCF:       $VCF"
    echo "Consensus: $CONS"

    bcftools consensus -f "$REF" "$VCF" > "$CONS"

    echo "Finished consensus for $SAMPLE"
done < "$SAMPLE_LIST"

echo "ALL CONSENSUS FASTAS FINISHED"
