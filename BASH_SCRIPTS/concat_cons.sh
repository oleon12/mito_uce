
set -euo pipefail

cd /scratch/odl7/sturnira_mito

KEEP_LIST="results/keep_samples_intersection_le40.txt"
OUTFILE="results/combined_masked_consensus_intersection.fasta"
FAILED_LIST="results/combined_masked_consensus_missing.txt"

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

# Ingroup
while read -r SAMPLE; do
    [[ -z "${SAMPLE:-}" ]] && continue
    FASTA="results/$SAMPLE/consensus_masked/${SAMPLE}_mt_masked.fasta"

    if [[ -f "$FASTA" ]]; then
        cat "$FASTA" >> "$OUTFILE"
        echo >> "$OUTFILE"
    else
        echo "$SAMPLE" >> "$FAILED_LIST"
    fi
done < "$KEEP_LIST"

# Outgroups
for SAMPLE in "${OUTGROUPS[@]}"; do
    FASTA="results/$SAMPLE/consensus_masked/${SAMPLE}_mt_masked.fasta"

    if [[ -f "$FASTA" ]]; then
        cat "$FASTA" >> "$OUTFILE"
        echo >> "$OUTFILE"
    else
        echo "$SAMPLE" >> "$FAILED_LIST"
    fi
done

echo "Combined masked consensus FASTA written to: $OUTFILE"
echo "Missing samples listed in: $FAILED_LIST"
echo "N sequences:" $(grep -c "^>" "$OUTFILE")
