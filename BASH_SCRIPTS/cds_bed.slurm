#!/bin/bash
#SBATCH --job-name=cds_bed
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=244G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira/LOG/cds_bed_%j.out
#SBATCH --error=/scratch/odl7/sturnira/ERR_OUT/cds_bed_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail

cd /scratch/odl7/sturnira_mito

############################################
# CONFIG
############################################
REF="references/S_ludovici_QCAZ_18312.fasta"
REFGB="references/S_ludovici_QCAZ_18312.gb"
OUTDIR="results"
OUTBED="$OUTDIR/cds_coords.bed"

############################################
# ABSOLUTE PATHS
############################################
REF=$(realpath "$REF")
REFGB=$(realpath "$REFGB")
OUTDIR=$(realpath -m "$OUTDIR")
OUTBED=$(realpath -m "$OUTBED")

############################################
# CHECK INPUTS
############################################
if [[ ! -f "$REF" ]]; then
    echo "ERROR: FASTA not found: $REF"
    exit 1
fi

if [[ ! -f "$REFGB" ]]; then
    echo "ERROR: GenBank file not found: $REFGB"
    exit 1
fi

mkdir -p "$OUTDIR"

############################################
# GET CHROM NAME FROM FASTA
############################################
CHROM=$(grep "^>" "$REF" | sed 's/>//' | awk '{print $1}')

if [[ -z "$CHROM" ]]; then
    echo "ERROR: Could not read sequence name from FASTA header."
    exit 1
fi

echo "Reference FASTA: $REF"
echo "Reference GB:    $REFGB"
echo "Chrom name:      $CHROM"
echo "Output BED:      $OUTBED"

############################################
# BUILD BED FROM CDS FEATURES
############################################
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
' "$REFGB" > "$OUTBED"

############################################
# CHECK OUTPUT
############################################
if [[ ! -s "$OUTBED" ]]; then
    echo "ERROR: BED file was created but is empty: $OUTBED"
    exit 1
fi

echo "CDS BED created successfully."
echo "Preview:"
head "$OUTBED"
