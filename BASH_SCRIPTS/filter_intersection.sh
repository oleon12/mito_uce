
set -euo pipefail

cd /scratch/odl7/sturnira_mito

CDS_KEEP="results/keep_samples_cds_le40.txt"
CONS_KEEP="results/keep_samples_masked_consensus_le40.txt"
OUT_INTERSECTION="results/keep_samples_intersection_le40.txt"

for f in "$CDS_KEEP" "$CONS_KEEP"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: file not found: $f"
        exit 1
    fi
done

sort -u "$CDS_KEEP" > /tmp/cds_keep_sorted_$$.txt
sort -u "$CONS_KEEP" > /tmp/cons_keep_sorted_$$.txt

comm -12 /tmp/cds_keep_sorted_$$.txt /tmp/cons_keep_sorted_$$.txt > "$OUT_INTERSECTION"

rm -f /tmp/cds_keep_sorted_$$.txt /tmp/cons_keep_sorted_$$.txt

echo "Intersection keep list written to: $OUT_INTERSECTION"
echo "Intersection count:" $(wc -l < "$OUT_INTERSECTION")
