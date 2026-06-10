
set -euo pipefail

cd /scratch/odl7/sturnira_mito

CONDA_ENV="mt_pipeline"
OUTGROUP_DIR="raw_data/outgroups"
RESULTS_DIR="results"
FAILED_LIST="$RESULTS_DIR/outgroup_extraction_failed.txt"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

OUTGROUP_DIR=$(realpath "$OUTGROUP_DIR")
RESULTS_DIR=$(realpath -m "$RESULTS_DIR")

: > "$FAILED_LIST"

for FASTA in "$OUTGROUP_DIR"/*.fasta; do
    [[ -f "$FASTA" ]] || continue

    BASENAME=$(basename "$FASTA" .fasta)
    GB="$OUTGROUP_DIR/${BASENAME}.gb"

    if [[ ! -f "$GB" ]]; then
        echo "Missing matching GB for: $BASENAME" | tee -a "$FAILED_LIST"
        continue
    fi

    OUTGROUP="$BASENAME"
    OUTDIR="$RESULTS_DIR/$OUTGROUP"
    CONSDIR="$OUTDIR/consensus_masked"
    GENEDIR="$OUTDIR/genes_masked"
    BED="$OUTDIR/${OUTGROUP}_cds_coords.bed"

    mkdir -p "$CONSDIR" "$GENEDIR"

    echo "========================================"
    echo "Processing outgroup: $OUTGROUP"

    cp "$FASTA" "$CONSDIR/${OUTGROUP}_mt_masked.fasta" || {
        echo "Copy FASTA failed: $OUTGROUP" | tee -a "$FAILED_LIST"
        continue
    }

    CHROM=$(grep "^>" "$FASTA" | sed 's/^>//' | awk '{print $1}')

    if [[ -z "$CHROM" ]]; then
        echo "Could not read FASTA header: $OUTGROUP" | tee -a "$FAILED_LIST"
        continue
    fi

    awk -v chrom="$CHROM" '
    /^     CDS/ {
        line=$0
        sub(/^     CDS[[:space:]]+/, "", line)
        strand="+"
        if (line ~ /complement/) strand="-"
        gsub(/complement\(|join\(|\)| /, "", line)
        n=split(line, parts, ",")
        for (i=1; i<=n; i++) {
            split(parts[i], coords, /\.\./)
            if (coords[1] ~ /^[0-9]+$/ && coords[2] ~ /^[0-9]+$/) {
                start=coords[1]-1
                end=coords[2]
                print chrom "\t" start "\t" end "\tCDS_part" i "\t0\t" strand
            }
        }
    }
    ' "$GB" > "$BED"

    if [[ ! -s "$BED" ]]; then
        echo "Empty CDS BED: $OUTGROUP" | tee -a "$FAILED_LIST"
        continue
    fi

    bedtools getfasta \
        -fi "$CONSDIR/${OUTGROUP}_mt_masked.fasta" \
        -bed "$BED" \
        -s \
        -fo "$GENEDIR/${OUTGROUP}_cds_masked.fasta" || {
        echo "bedtools getfasta failed: $OUTGROUP" | tee -a "$FAILED_LIST"
        continue
    }

    TMP="$GENEDIR/${OUTGROUP}_cds_masked.fasta.tmp"
    awk -v sample="$OUTGROUP" '
        /^>/ {print ">" sample "|" ++i; next}
        {print}
    ' "$GENEDIR/${OUTGROUP}_cds_masked.fasta" > "$TMP"
    mv "$TMP" "$GENEDIR/${OUTGROUP}_cds_masked.fasta"

    TMP2="$CONSDIR/${OUTGROUP}_mt_masked.fasta.tmp"
    awk -v sample="$OUTGROUP" '
        BEGIN {done=0}
        /^>/ {
            if (done==0) {
                print ">" sample
                done=1
            } else {
                print
            }
            next
        }
        {print}
    ' "$CONSDIR/${OUTGROUP}_mt_masked.fasta" > "$TMP2"
    mv "$TMP2" "$CONSDIR/${OUTGROUP}_mt_masked.fasta"

    echo "Finished: $OUTGROUP"
done

echo "Done. Failures listed in: $FAILED_LIST"
