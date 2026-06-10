#!/bin/bash
#SBATCH --job-name=vcf_mt_all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/vcf_mt_all_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/vcf_mt_all_%j.err
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

    MAPPED_BAM="$OUTDIR/$SAMPLE/bam/$SAMPLE.mapped.bam"
    VCFDIR="$OUTDIR/$SAMPLE/vcf"
    VCF="$VCFDIR/$SAMPLE.vcf.gz"

    if [[ ! -f "$MAPPED_BAM" ]]; then
        echo "WARNING: mapped BAM not found for $SAMPLE -> $MAPPED_BAM"
        continue
    fi

    mkdir -p "$VCFDIR"

    if [[ -f "$VCF" && ( -f "${VCF}.csi" || -f "${VCF}.tbi" ) ]]; then
        echo "VCF already exists for $SAMPLE. Skipping."
        continue
    fi

    echo "========================================"
    echo "Processing sample: $SAMPLE"
    echo "Mapped BAM: $MAPPED_BAM"
    echo "VCF:        $VCF"

    bcftools mpileup -f "$REF" "$MAPPED_BAM" -Ou | \
        bcftools call -mv --ploidy 1 -Oz -o "$VCF"

    bcftools index "$VCF"

    echo "Finished variant calling for $SAMPLE"
done < "$SAMPLE_LIST"

echo "ALL SAMPLES FINISHED"
