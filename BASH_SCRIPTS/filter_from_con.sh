

set -euo pipefail

cd /scratch/odl7/sturnira_mito

INPUT="results/masked_consensus_stats.tsv"
OUT_KEEP="results/keep_samples_masked_consensus_le40.txt"
OUT_DROP="results/drop_samples_masked_consensus_gt40.txt"
CUTOFF=40

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT"
    exit 1
fi

awk -F'\t' -v cutoff="$CUTOFF" '
NR==1 {next}
{
    sample=$1
    pct=$5 + 0
    if (pct <= cutoff) print sample
}
' "$INPUT" > "$OUT_KEEP"

awk -F'\t' -v cutoff="$CUTOFF" '
NR==1 {next}
{
    sample=$1
    pct=$5 + 0
    if (pct > cutoff) print sample
}
' "$INPUT" > "$OUT_DROP"

echo "Masked consensus keep list written to: $OUT_KEEP"
echo "Masked consensus drop list written to: $OUT_DROP"
echo "Kept:" $(wc -l < "$OUT_KEEP")
echo "Dropped:" $(wc -l < "$OUT_DROP")
