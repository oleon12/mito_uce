#!/usr/bin/env bash
#SBATCH --job-name=cds_bed
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/cds_bed_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/cds_bed_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

REF="references/S_ludovici_QCAZ_18312.fasta"
REFGB="references/S_ludovici_QCAZ_18312.gb"

OUTDIR="results"
# Coding portion only: terminal complete/incomplete stop bases removed.
OUTBED="$OUTDIR/cds_coords.bed"
# Full CDS coordinates exactly as annotated in GenBank.
FULLBED="$OUTDIR/cds_coords_annotated_full.bed"
# Auditable annotation and validation report.
OUTTSV="$OUTDIR/cds_coords.tsv"
# Genomic overlaps among retained coding regions; important before concatenation.
OUTOVERLAP="$OUTDIR/cds_overlaps.tsv"

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found in Conda environment: $CONDA_ENV" >&2
    exit 1
}

python3 - <<'PY_CHECK'
try:
    import Bio
except ImportError as exc:
    raise SystemExit(
        "ERROR: Biopython is required in the active environment. "
        "Install it with: conda install -c conda-forge biopython"
    ) from exc
PY_CHECK

############################################
# RESOLVE AND CHECK PATHS
############################################
REF=$(realpath "$REF")
REFGB=$(realpath "$REFGB")
OUTDIR=$(realpath -m "$OUTDIR")
OUTBED=$(realpath -m "$OUTBED")
FULLBED=$(realpath -m "$FULLBED")
OUTTSV=$(realpath -m "$OUTTSV")
OUTOVERLAP=$(realpath -m "$OUTOVERLAP")

[[ -s "$REF" ]] || {
    echo "ERROR: Reference FASTA not found or empty: $REF" >&2
    exit 1
}

[[ -s "$REFGB" ]] || {
    echo "ERROR: Reference GenBank not found or empty: $REFGB" >&2
    exit 1
}

mkdir -p "$OUTDIR"

############################################
# PARSE AND VALIDATE THE GENBANK ANNOTATION
############################################
python3 - "$REF" "$REFGB" "$OUTBED" "$FULLBED" "$OUTTSV" "$OUTOVERLAP" <<'PY'
import csv
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

ref_fasta = Path(sys.argv[1])
ref_genbank = Path(sys.argv[2])
coding_bed = Path(sys.argv[3])
full_bed = Path(sys.argv[4])
metadata_tsv = Path(sys.argv[5])
overlap_tsv = Path(sys.argv[6])

EXPECTED_GENES = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"
]

GENE_ALIASES = {
    "COI": "COX1",
    "CO1": "COX1",
    "COII": "COX2",
    "CO2": "COX2",
    "COIII": "COX3",
    "CO3": "COX3",
    "CYB": "CYTB",
    "COB": "CYTB",
    "ATPASE6": "ATP6",
    "ATPASE8": "ATP8",
    "NAD1": "ND1",
    "NAD2": "ND2",
    "NAD3": "ND3",
    "NAD4": "ND4",
    "NAD4L": "ND4L",
    "NAD5": "ND5",
    "NAD6": "ND6",
}


def canonical_gene_name(value):
    normalized = (
        value.upper()
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
    )
    if normalized in EXPECTED_GENES:
        return normalized
    return GENE_ALIASES.get(normalized, normalized)


def one_record(path, fmt):
    records = list(SeqIO.parse(str(path), fmt))
    if len(records) != 1:
        raise SystemExit(
            "ERROR: {} must contain exactly one sequence record; found {}."
            .format(path, len(records))
        )
    return records[0]


def genomic_blocks(location):
    """Return zero-based, half-open blocks sorted by genomic coordinate."""
    parts = list(getattr(location, "parts", [location]))
    blocks = sorted(
        [(int(part.start), int(part.end)) for part in parts],
        key=lambda x: (x[0], x[1])
    )
    if not blocks or any(start >= end for start, end in blocks):
        raise SystemExit("ERROR: invalid or empty CDS location: {}".format(location))
    return blocks


def retain_first_coding_bases(blocks, strand, coding_length):
    """
    Retain the first coding_length nucleotides in transcript orientation.

    For + features, transcript order is low to high.
    For - features, transcript order is high to low.
    Returned blocks are sorted in ascending genomic order for BED12.
    """
    ordered = blocks if strand == "+" else list(reversed(blocks))
    remaining = coding_length
    kept = []

    for start, end in ordered:
        if remaining <= 0:
            break

        block_length = end - start
        take = min(block_length, remaining)

        if strand == "+":
            kept.append((start, start + take))
        else:
            kept.append((end - take, end))

        remaining -= take

    if remaining != 0:
        raise SystemExit(
            "ERROR: CDS blocks contain fewer bases than the requested "
            "coding length ({} bases missing).".format(remaining)
        )

    return sorted(kept, key=lambda x: (x[0], x[1]))


def bed12_line(chrom, gene, strand, blocks):
    chrom_start = min(start for start, _ in blocks)
    chrom_end = max(end for _, end in blocks)
    block_sizes = [end - start for start, end in blocks]
    block_starts = [start - chrom_start for start, _ in blocks]

    fields = [
        chrom,
        str(chrom_start),
        str(chrom_end),
        gene,
        "0",
        strand,
        str(chrom_start),
        str(chrom_end),
        "0",
        str(len(blocks)),
        ",".join(map(str, block_sizes)) + ",",
        ",".join(map(str, block_starts)) + ",",
    ]
    return "\t".join(fields)


fasta_record = one_record(ref_fasta, "fasta")
gb_record = one_record(ref_genbank, "genbank")

if len(fasta_record.seq) != len(gb_record.seq):
    raise SystemExit(
        "ERROR: FASTA length ({}) differs from GenBank length ({})."
        .format(len(fasta_record.seq), len(gb_record.seq))
    )

if str(fasta_record.seq).upper() != str(gb_record.seq).upper():
    raise SystemExit(
        "ERROR: FASTA and GenBank nucleotide sequences are not identical. "
        "GenBank coordinates cannot safely be applied to this FASTA."
    )

chrom = fasta_record.id
features = {}

for feature in gb_record.features:
    if feature.type != "CDS":
        continue

    raw_name = (
        feature.qualifiers.get("gene")
        or feature.qualifiers.get("locus_tag")
        or feature.qualifiers.get("product")
        or [""]
    )[0]
    gene = canonical_gene_name(raw_name)

    if gene not in EXPECTED_GENES:
        continue

    if gene in features:
        raise SystemExit(
            "ERROR: duplicate CDS annotation for gene {}.".format(gene)
        )

    features[gene] = feature

missing = [gene for gene in EXPECTED_GENES if gene not in features]
extra = sorted(set(features) - set(EXPECTED_GENES))

if missing:
    raise SystemExit(
        "ERROR: missing expected mitochondrial CDS annotations: {}."
        .format(", ".join(missing))
    )
if extra:
    raise SystemExit(
        "ERROR: unexpected CDS annotations were retained: {}."
        .format(", ".join(extra))
    )

metadata_rows = []
coding_lines = []
full_lines = []
coding_blocks_by_gene = {}

for gene in EXPECTED_GENES:
    feature = features[gene]

    strand_value = feature.location.strand
    if strand_value not in (1, -1):
        raise SystemExit(
            "ERROR: gene {} has undefined strand: {}."
            .format(gene, strand_value)
        )
    strand = "+" if strand_value == 1 else "-"

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    if codon_start != 1:
        raise SystemExit(
            "ERROR: gene {} has codon_start={}. This workflow expects a "
            "complete mitochondrial reference with codon_start=1."
            .format(gene, codon_start)
        )

    transl_table = int(feature.qualifiers.get("transl_table", ["2"])[0])
    translation = (
        feature.qualifiers.get("translation", [""])[0]
        .replace(" ", "")
        .replace("\n", "")
    )
    if not translation:
        raise SystemExit(
            "ERROR: gene {} lacks an annotated /translation qualifier."
            .format(gene)
        )

    full_blocks = genomic_blocks(feature.location)
    annotated_nt = sum(end - start for start, end in full_blocks)
    coding_nt = len(translation) * 3
    terminal_stop_nt = annotated_nt - coding_nt

    if terminal_stop_nt not in (0, 1, 2, 3):
        raise SystemExit(
            "ERROR: gene {} has annotated length {} but translation implies "
            "{} coding nucleotides; difference={} is not a valid terminal "
            "stop-codon length."
            .format(gene, annotated_nt, coding_nt, terminal_stop_nt)
        )

    coding_blocks = retain_first_coding_bases(
        full_blocks, strand, coding_nt
    )

    # Independent sequence-level validation using Biopython's feature parser.
    full_sequence = feature.extract(gb_record.seq)
    coding_sequence = full_sequence[:coding_nt]
    translated = str(coding_sequence.translate(table=transl_table))

    # GenBank CDS translations represent a valid alternative initiation codon
    # as methionine. Plain Seq.translate() reports the codon's internal amino
    # acid, so normalize the first residue when the codon table recognizes it
    # as an allowed start codon.
    codon_table = CodonTable.unambiguous_dna_by_id[transl_table]
    first_codon = str(coding_sequence[:3]).upper()
    if first_codon in codon_table.start_codons and translated:
        translated = "M" + translated[1:]

    if translated != translation:
        raise SystemExit(
            "ERROR: translated coding sequence for {} does not match the "
            "GenBank /translation annotation.".format(gene)
        )

    product = feature.qualifiers.get("product", [""])[0]
    start_1based = min(start for start, _ in full_blocks) + 1
    end_1based = max(end for _, end in full_blocks)

    full_lines.append(bed12_line(chrom, gene, strand, full_blocks))
    coding_lines.append(bed12_line(chrom, gene, strand, coding_blocks))
    coding_blocks_by_gene[gene] = coding_blocks

    metadata_rows.append({
        "gene": gene,
        "product": product,
        "strand": strand,
        "start_1based": start_1based,
        "end_1based": end_1based,
        "block_count": len(full_blocks),
        "annotated_nt": annotated_nt,
        "coding_nt_retained": coding_nt,
        "terminal_stop_nt_removed": terminal_stop_nt,
        "protein_aa": len(translation),
        "codon_start": codon_start,
        "transl_table": transl_table,
        "translation_validated": "yes",
    })

full_bed.write_text("\n".join(full_lines) + "\n")
coding_bed.write_text("\n".join(coding_lines) + "\n")

fieldnames = [
    "gene",
    "product",
    "strand",
    "start_1based",
    "end_1based",
    "block_count",
    "annotated_nt",
    "coding_nt_retained",
    "terminal_stop_nt_removed",
    "protein_aa",
    "codon_start",
    "transl_table",
    "translation_validated",
]

with metadata_tsv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(metadata_rows)

# Record genomic overlaps so duplicated nucleotide positions are handled
# explicitly when the gene alignments are concatenated later.
overlap_rows = []
for i, gene1 in enumerate(EXPECTED_GENES):
    for gene2 in EXPECTED_GENES[i + 1:]:
        overlap_nt = 0
        for start1, end1 in coding_blocks_by_gene[gene1]:
            for start2, end2 in coding_blocks_by_gene[gene2]:
                overlap_nt += max(0, min(end1, end2) - max(start1, start2))
        if overlap_nt > 0:
            overlap_rows.append((gene1, gene2, overlap_nt))

with overlap_tsv.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene1", "gene2", "coding_overlap_nt"])
    writer.writerows(overlap_rows)

print("Reference FASTA record: {}".format(fasta_record.id))
print("GenBank record:         {}".format(gb_record.id))
print("Reference length:        {} bp".format(len(fasta_record.seq)))
print("Validated CDS genes:     {}".format(len(metadata_rows)))
print("Retained coding length:  {} bp".format(
    sum(row["coding_nt_retained"] for row in metadata_rows)
))
print("Coding BED12:            {}".format(coding_bed))
print("Full annotation BED12:   {}".format(full_bed))
print("Metadata/QC table:       {}".format(metadata_tsv))
print("Coding-overlap report:   {}".format(overlap_tsv))
PY

############################################
# FINAL OUTPUT CHECKS
############################################
for file in "$OUTBED" "$FULLBED" "$OUTTSV" "$OUTOVERLAP"; do
    [[ -s "$file" ]] || {
        echo "ERROR: expected output is missing or empty: $file" >&2
        exit 1
    }
done

BED_GENES=$(cut -f4 "$OUTBED" | wc -l)
UNIQUE_GENES=$(cut -f4 "$OUTBED" | sort -u | wc -l)

if [[ "$BED_GENES" -ne 13 || "$UNIQUE_GENES" -ne 13 ]]; then
    echo "ERROR: expected 13 unique genes in $OUTBED; found $BED_GENES rows and $UNIQUE_GENES unique names." >&2
    exit 1
fi

echo
echo "CDS annotation completed successfully."
echo "The primary BED is BED12 and contains coding sequence only."
echo "Later extraction must use: bedtools getfasta -split -s -name"
echo
column -t -s $'\t' "$OUTTSV"
echo
echo "Retained coding-region overlaps:"
column -t -s $'\t' "$OUTOVERLAP"
