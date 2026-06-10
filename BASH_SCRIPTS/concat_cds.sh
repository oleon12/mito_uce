
set -euo pipefail

cd /scratch/odl7/sturnira_mito

KEEP_LIST="results/keep_samples_intersection_le40.txt"
OUTFILE="results/combined_masked_cds_intersection.fasta"
FAILED_LIST="results/combined_masked_cds_missing.txt"

# Add outgroups here, space-separated
OUTGROUPS=(
  "Artibeus_PP853570.1"
  "Glossophaga_NC_065682.1"
)

: > "$OUTFILE"
: > "$FAILED_LIST"

if [[ ! -f "$KEEP_LIST" ]]; then
    echo "ERROR: keep list not found: $KEEP_LIST"
    exit 1
fi

write_one_sample() {
    local sample="$1"
    local fasta="$2"

    if [[ ! -f "$fasta" ]]; then
        echo "$sample" >> "$FAILED_LIST"
        return
    fi

    seq=$(grep -v "^>" "$fasta" | tr -d '\n\r')

    if [[ -z "$seq" ]]; then
        echo "$sample" >> "$FAILED_LIST"
        return
    fi

    {
        echo ">$sample"
        echo "$seq"
    } >> "$OUTFILE"
}

# Ingroup
while read -r SAMPLE; do
    [[ -z "${SAMPLE:-}" ]] && continue
    FASTA="results/$SAMPLE/genes_masked/${SAMPLE}_cds_masked.fasta"
    write_one_sample "$SAMPLE" "$FASTA"
done < "$KEEP_LIST"

# Outgroups
for SAMPLE in "${OUTGROUPS[@]}"; do
    FASTA="results/$SAMPLE/genes_masked/${SAMPLE}_cds_masked.fasta"
    write_one_sample "$SAMPLE" "$FASTA"
done

echo "Combined masked CDS FASTA written to: $OUTFILE"
echo "Missing samples listed in: $FAILED_LIST"
echo "N sequences:" $(grep -c "^>" "$OUTFILE")
