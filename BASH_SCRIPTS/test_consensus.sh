set -eo pipefail

cd /scratch/odl7/sturnira_mito

OUTFILE="results/consensus_stats.tsv"

echo -e "Sample\tFASTA\tLength\tNs\tPercent_N" > "$OUTFILE"

for fasta in results/*/consensus/*_mt.fasta; do
    [[ -f "$fasta" ]] || continue

    sample=$(basename "$fasta" _mt.fasta)

    seq=$(grep -v "^>" "$fasta" | tr -d '\n\r')
    len=${#seq}

    ncount=$(printf "%s" "$seq" | grep -o "N" | wc -l || true)

    if [[ "$len" -gt 0 ]]; then
        pct=$(awk -v n="$ncount" -v l="$len" 'BEGIN {printf "%.4f", (n/l)*100}')
    else
        pct="0.0000"
    fi

    echo -e "${sample}\t${fasta}\t${len}\t${ncount}\t${pct}" >> "$OUTFILE"
done

echo "Consensus stats written to: $OUTFILE"
