

set -euo pipefail

cd /scratch/odl7/sturnira_mito

SAMPLE_LIST="CONFS/sample_list.txt"
OUTDIR="results"
FAILED_LIST="$OUTDIR/rename_fasta_headers_failed.txt"

: > "$FAILED_LIST"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: sample list not found: $SAMPLE_LIST"
    exit 1
fi

while read -r SAMPLE READ1 READ2; do
    [[ -z "${SAMPLE:-}" ]] && continue

    echo "========================================"
    echo "Processing: $SAMPLE"

    ############################################
    # 1) masked consensus fasta
    ############################################
    CONS="$OUTDIR/$SAMPLE/consensus_masked/${SAMPLE}_mt_masked.fasta"

    if [[ -f "$CONS" ]]; then
        TMP="${CONS}.tmp"
        awk -v sample="$SAMPLE" '
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
        ' "$CONS" > "$TMP" && mv "$TMP" "$CONS" || {
            echo "Failed renaming header in masked consensus: $SAMPLE" | tee -a "$FAILED_LIST"
            rm -f "$TMP"
        }
    else
        echo "Missing masked consensus fasta: $SAMPLE" | tee -a "$FAILED_LIST"
    fi

    ############################################
    # 2) masked CDS fasta
    ############################################
    CDS="$OUTDIR/$SAMPLE/genes_masked/${SAMPLE}_cds_masked.fasta"

    if [[ -f "$CDS" ]]; then
        TMP="${CDS}.tmp"
        awk -v sample="$SAMPLE" '
            /^>/ {print ">" sample "|" ++i; next}
            {print}
        ' "$CDS" > "$TMP" && mv "$TMP" "$CDS" || {
            echo "Failed renaming header in masked CDS: $SAMPLE" | tee -a "$FAILED_LIST"
            rm -f "$TMP"
        }
    else
        echo "Missing masked CDS fasta: $SAMPLE" | tee -a "$FAILED_LIST"
    fi

done < "$SAMPLE_LIST"

echo "Done. Failures listed in: $FAILED_LIST"
