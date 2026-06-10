#!/bin/bash
#SBATCH --job-name=masked_cons_all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/masked_cons_all_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/masked_cons_all_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/odl7/sturnira_mito

############################################
# CONFIG
############################################
CONDA_ENV="mt_pipeline"
SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
OUTDIR="results"
MIN_DEPTH=3

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

FAILED_LIST="$OUTDIR/masked_consensus_failed_samples.txt"
: > "$FAILED_LIST"

while read -r SAMPLE READ1 READ2; do
    [[ -z "$SAMPLE" ]] && continue

    BAM="$OUTDIR/$SAMPLE/bam/$SAMPLE.mapped.bam"
    VCF="$OUTDIR/$SAMPLE/vcf/$SAMPLE.vcf.gz"
    CONSDIR="$OUTDIR/$SAMPLE/consensus_masked"
    LOWBED="$CONSDIR/${SAMPLE}.lowcov.bed"
    CLIPPEDBED="$CONSDIR/${SAMPLE}.lowcov.clipped.bed"
    RAWCONS="$CONSDIR/${SAMPLE}.raw_consensus.fa"
    CONS="$CONSDIR/${SAMPLE}_mt_masked.fasta"
    TMPCONS="$CONSDIR/${SAMPLE}_mt_masked.fasta.tmp"

    mkdir -p "$CONSDIR"
    rm -f "$TMPCONS"

    if [[ -f "$CONS" ]]; then
        echo "Masked consensus already exists for $SAMPLE. Skipping."
        continue
    fi

    echo "========================================"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$BAM" ]]; then
        echo "Missing BAM: $SAMPLE" | tee -a "$FAILED_LIST"
        continue
    fi

    if [[ ! -f "$VCF" ]]; then
        echo "Missing VCF: $SAMPLE" | tee -a "$FAILED_LIST"
        continue
    fi

    ############################################
    # 1) BUILD RAW CONSENSUS
    ############################################
    nvar=$(bcftools view -H "$VCF" | wc -l)

    if [[ "$nvar" -eq 0 ]]; then
        echo "No variants in VCF for $SAMPLE; using reference as raw consensus"
        cp "$REF" "$RAWCONS" || {
            echo "Copy reference failed: $SAMPLE" | tee -a "$FAILED_LIST"
            continue
        }
    else
        bcftools consensus -f "$REF" "$VCF" > "$RAWCONS" || {
            echo "bcftools consensus failed: $SAMPLE" | tee -a "$FAILED_LIST"
            rm -f "$RAWCONS"
            continue
        }
    fi

    ############################################
    # 2) BUILD LOW-COVERAGE BED FROM BAM
    ############################################
    samtools depth -aa "$BAM" | awk -v min="$MIN_DEPTH" '
        $3 < min {print $1 "\t" ($2-1) "\t" $2}
    ' > "$LOWBED" || {
        echo "Depth failed: $SAMPLE" | tee -a "$FAILED_LIST"
        rm -f "$RAWCONS"
        continue
    }

    ############################################
    # 3) CLIP BED TO RAW CONSENSUS LENGTH
    ############################################
    rawlen=$(grep -v "^>" "$RAWCONS" | tr -d '\n\r' | wc -c)

    awk -v max="$rawlen" 'BEGIN{OFS="\t"}
        $2 < max {
            if ($3 > max) $3=max
            if ($3 > $2) print $1,$2,$3
        }
    ' "$LOWBED" > "$CLIPPEDBED" || {
        echo "BED clipping failed: $SAMPLE" | tee -a "$FAILED_LIST"
        rm -f "$RAWCONS"
        continue
    }

    ############################################
    # 4) IF NO MASKING NEEDED, COPY RAW CONSENSUS
    ############################################
    if [[ ! -s "$CLIPPEDBED" ]]; then
        cp "$RAWCONS" "$CONS" || {
            echo "Copy raw consensus failed: $SAMPLE" | tee -a "$FAILED_LIST"
            rm -f "$RAWCONS"
            continue
        }
        echo "No low-coverage sites after clipping for $SAMPLE; copied raw consensus"
        continue
    fi

    ############################################
    # 5) MASK LOW-COVERAGE POSITIONS
    ############################################
    bedtools maskfasta \
        -fi "$RAWCONS" \
        -bed "$CLIPPEDBED" \
        -fo "$TMPCONS" || {
        echo "maskfasta failed: $SAMPLE" | tee -a "$FAILED_LIST"
        rm -f "$TMPCONS"
        continue
    }

    ############################################
    # 6) CHECK LENGTH AND FINALIZE
    ############################################
    seqlen=$(grep -v "^>" "$TMPCONS" | tr -d '\n\r' | wc -c)

    if [[ "$seqlen" -lt 15000 ]]; then
        echo "Short consensus (<15000 bp): $SAMPLE length=$seqlen" | tee -a "$FAILED_LIST"
        rm -f "$TMPCONS"
        continue
    fi

    mv "$TMPCONS" "$CONS"
    echo "Finished: $SAMPLE length=$seqlen"

done < "$SAMPLE_LIST"

echo "Done. Failed samples listed in: $FAILED_LIST"
