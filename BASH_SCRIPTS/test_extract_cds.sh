
set -uo pipefail

cd /scratch/odl7/sturnira_mito

OUTFILE="results/masked_cds_summary.tsv"

echo -e "Sample\tFASTA\tN_CDS\tTotal_Length\tNs\tPercent_N" > "$OUTFILE"

for fasta in results/*/genes_masked/*_cds_masked.fasta; do
    [[ -f "$fasta" ]] || continue

    sample=$(basename "$fasta" _cds_masked.fasta)
    nseq=$(grep -c "^>" "$fasta" 2>/dev/null || echo 0)

    seq=$(grep -v "^>" "$fasta" | tr -d '\n\r')
    len=${#seq}
    ncount=$(printf "%s" "$seq" | grep -o "N" | wc -l || true)

    if [[ "$len" -gt 0 ]]; then
        pct=$(awk -v n="$ncount" -v l="$len" 'BEGIN {printf "%.4f", (n/l)*100}')
    else
        pct="0.0000"
    fi

    echo -e "${sample}\t${fasta}\t${nseq}\t${len}\t${ncount}\t${pct}" >> "$OUTFILE"
done

echo "Masked CDS summary written to: $OUTFILE"
