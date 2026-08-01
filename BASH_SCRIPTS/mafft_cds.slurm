#!/usr/bin/env bash
#SBATCH --job-name=mafft_13pcg
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/mafft_13pcg_%A_%a.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/mafft_13pcg_%A_%a.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Align each of the 13 mitochondrial protein-coding genes independently:
#
# nucleotide CDS
#   -> vertebrate-mitochondrial protein translation
#   -> MAFFT protein alignment
#   -> codon-aware back-translation
#   -> concatenated nucleotide supermatrix
#   -> 39 gene-by-codon-position partitions
#
# This script processes M50, M40, M30, and M20 as a SLURM array. It consumes the
# unaligned gene matrices produced by the revised concat_cds.sh.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

THRESHOLDS=(50 40 30 20)
THRESHOLD="${THRESHOLDS[$SLURM_ARRAY_TASK_ID]}"
LABEL="M${THRESHOLD}"

INPUT_ROOT="results/cds_gene_matrices"
INPUT_DIR="$INPUT_ROOT/$LABEL/unaligned_by_gene"
TAXA_FILE="$INPUT_ROOT/$LABEL/taxa.txt"

# Coding-only BED12 produced by the revised cds_bed.slurm. It is used here
# to document genomic overlaps that are represented twice when full genes
# are concatenated.
CDS_BED="results/cds_coords.bed"

# M50 determines one common MAFFT strategy for all four sensitivity datasets.
M50_TAXA_FILE="$INPUT_ROOT/M50/taxa.txt"

OUTROOT="results/alignments/cds"
OUTDIR="$OUTROOT/$LABEL"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

EXPECTED_GENES=(
    "ND1"
    "ND2"
    "COX1"
    "COX2"
    "ATP8"
    "ATP6"
    "COX3"
    "ND3"
    "ND4L"
    "ND4"
    "ND5"
    "ND6"
    "CYTB"
)

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 mafft; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: $command_name was not found in $CONDA_ENV." >&2
        exit 1
    }
done

python3 - <<'PY_CHECK'
import sys

if sys.version_info < (3, 6):
    raise SystemExit("ERROR: Python 3.6 or newer is required.")

try:
    import Bio
except ImportError as exc:
    raise SystemExit(
        "ERROR: Biopython is required. Install it with: "
        "conda install -c conda-forge biopython"
    ) from exc
PY_CHECK

echo "MAFFT version:"
mafft --version 2>&1

############################################
# RESOLVE AND CHECK PATHS
############################################
INPUT_DIR=$(realpath "$INPUT_DIR")
TAXA_FILE=$(realpath "$TAXA_FILE")
M50_TAXA_FILE=$(realpath "$M50_TAXA_FILE")
CDS_BED=$(realpath "$CDS_BED")
OUTROOT=$(realpath -m "$OUTROOT")
OUTDIR=$(realpath -m "$OUTDIR")

[[ -d "$INPUT_DIR" ]] || {
    echo "ERROR: input gene directory not found: $INPUT_DIR" >&2
    exit 1
}

for file in "$TAXA_FILE" "$M50_TAXA_FILE" "$CDS_BED"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required taxa list missing or empty: $file" >&2
        exit 1
    }
done

for gene in "${EXPECTED_GENES[@]}"; do
    [[ -s "$INPUT_DIR/${gene}.fasta" ]] || {
        echo "ERROR: missing gene matrix: $INPUT_DIR/${gene}.fasta" >&2
        exit 1
    }
done

mkdir -p "$OUTROOT"

GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

############################################
# SELECT ONE STRATEGY FOR ALL THRESHOLDS
############################################
M50_N_TAXA=$(grep -cv '^[[:space:]]*$' "$M50_TAXA_FILE" || true)

if [[ "$M50_N_TAXA" -le 200 ]]; then
    STRATEGY="G-INS-i"
else
    STRATEGY="FFT-NS-i-1000"
fi

echo "Dataset:             $LABEL"
echo "M50 total taxa:      $M50_N_TAXA"
echo "Common strategy:     $STRATEGY"
echo "Threads:             $THREADS"

############################################
# ALIGN, BACK-TRANSLATE, AND CONCATENATE
############################################
python3 - \
    "$INPUT_DIR" \
    "$TAXA_FILE" \
    "$CDS_BED" \
    "$OUTDIR" \
    "$THREADS" \
    "$STRATEGY" \
    "$GENES_TEXT" <<'PY'
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


input_dir = Path(sys.argv[1])
taxa_file = Path(sys.argv[2])
cds_bed = Path(sys.argv[3])
outdir = Path(sys.argv[4])
threads = int(sys.argv[5])
strategy = sys.argv[6]
gene_order = [x for x in sys.argv[7].split(",") if x]

allowed_nt = set("ACGTRYSWKMBDHVN")
unambiguous_nt = set("ACGT")
allowed_aa = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-?")


def read_taxa(path):
    taxa = []
    seen = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            taxon = line.split()[0]

            if taxon in seen:
                raise ValueError(
                    "Duplicate taxon '{}' in {} at line {}."
                    .format(taxon, path, line_number)
                )

            seen.add(taxon)
            taxa.append(taxon)

    if not taxa:
        raise ValueError("No taxa found in {}.".format(path))

    return taxa


def coding_overlap_report(path, expected_gene_order):
    """Return positive pairwise genomic overlaps from a coding-only BED12."""
    blocks_by_gene = {}
    observed_order = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 12:
                raise ValueError(
                    "{} line {} is not BED12."
                    .format(path, line_number)
                )

            chrom_start = int(fields[1])
            gene = fields[3].upper()
            block_count = int(fields[9])
            sizes = [
                int(value)
                for value in fields[10].rstrip(",").split(",")
                if value
            ]
            starts = [
                int(value)
                for value in fields[11].rstrip(",").split(",")
                if value
            ]

            if gene in blocks_by_gene:
                raise ValueError(
                    "Duplicate gene '{}' in {}."
                    .format(gene, path)
                )

            if block_count != len(sizes) or block_count != len(starts):
                raise ValueError(
                    "Inconsistent BED12 blocks for '{}' in {}."
                    .format(gene, path)
                )

            blocks = [
                (chrom_start + relative, chrom_start + relative + size)
                for size, relative in zip(sizes, starts)
            ]

            blocks_by_gene[gene] = blocks
            observed_order.append(gene)

    if observed_order != expected_gene_order:
        raise ValueError(
            "Gene order/content in {} does not match the alignment workflow."
            .format(path)
        )

    positions_by_gene = {}

    for gene, blocks in blocks_by_gene.items():
        positions = set()

        for start, end in blocks:
            positions.update(range(start, end))

        positions_by_gene[gene] = positions

    rows = []

    for index, gene1 in enumerate(expected_gene_order):
        for gene2 in expected_gene_order[index + 1:]:
            overlap = positions_by_gene[gene1] & positions_by_gene[gene2]

            if overlap:
                rows.append({
                    "gene1": gene1,
                    "gene2": gene2,
                    "overlap_nt": len(overlap),
                    "first_reference_position_1based": min(overlap) + 1,
                    "last_reference_position_1based": max(overlap) + 1,
                })

    return rows


def read_gene_matrix(path, expected_taxa):
    records = list(SeqIO.parse(str(path), "fasta"))

    if not records:
        raise ValueError("No FASTA records found in {}.".format(path))

    ids = [record.id for record in records]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate FASTA identifiers in {}.".format(path))

    if ids != expected_taxa:
        raise ValueError(
            "Taxon names/order in {} do not match {}."
            .format(path, taxa_file)
        )

    sequences = {}

    for record in records:
        sequence = str(record.seq).upper().replace("U", "T")

        if not sequence:
            raise ValueError(
                "Empty sequence for '{}' in {}."
                .format(record.id, path)
            )

        invalid = sorted(set(sequence) - allowed_nt)

        if invalid:
            raise ValueError(
                "Invalid nucleotide symbols for '{}' in {}: {}."
                .format(record.id, path, ",".join(invalid))
            )

        if "-" in sequence:
            raise ValueError(
                "{} already contains gaps. Input must be unaligned."
                .format(path)
            )

        if len(sequence) % 3 != 0:
            raise ValueError(
                "{} in {} has length {}, not divisible by 3."
                .format(record.id, path, len(sequence))
            )

        sequences[record.id] = sequence

    return sequences


def definite_stop_codons(sequence, table=2):
    positions = []

    for start in range(0, len(sequence), 3):
        codon = sequence[start:start + 3]

        if set(codon).issubset(unambiguous_nt):
            aa = str(Seq(codon).translate(table=table))

            if aa == "*":
                positions.append(start // 3 + 1)

    return positions


def protein_for_alignment(sequence):
    protein = str(Seq(sequence).translate(table=2))

    invalid = sorted(set(protein.upper()) - allowed_aa)

    if invalid:
        raise ValueError(
            "Unexpected translated amino-acid symbols: {}."
            .format(",".join(invalid))
        )

    return protein.upper()


def run_mafft(protein_input, protein_output, log_path):
    if strategy == "G-INS-i":
        command = [
            "mafft",
            "--globalpair",
            "--maxiterate",
            "1000",
            "--thread",
            str(threads),
            str(protein_input),
        ]
    elif strategy == "FFT-NS-i-1000":
        command = [
            "mafft",
            "--retree",
            "2",
            "--maxiterate",
            "1000",
            "--thread",
            str(threads),
            str(protein_input),
        ]
    else:
        raise ValueError("Unknown MAFFT strategy: {}".format(strategy))

    with protein_output.open("w") as stdout, log_path.open("w") as stderr:
        completed = subprocess.run(
            command,
            stdout=stdout,
            stderr=stderr,
        )

    if completed.returncode != 0:
        raise RuntimeError(
            "MAFFT failed for {}. See {}."
            .format(protein_input, log_path)
        )

    return " ".join(command)


def backtranslate(aligned_protein, nucleotide_sequence, taxon, gene):
    codons = [
        nucleotide_sequence[start:start + 3]
        for start in range(0, len(nucleotide_sequence), 3)
    ]

    aligned_codons = []
    codon_index = 0

    for amino_acid in aligned_protein:
        if amino_acid == "-":
            aligned_codons.append("---")
            continue

        if codon_index >= len(codons):
            raise ValueError(
                "Protein alignment consumes too many codons for {}|{}."
                .format(taxon, gene)
            )

        aligned_codons.append(codons[codon_index])
        codon_index += 1

    if codon_index != len(codons):
        raise ValueError(
            "Protein alignment consumed {} of {} codons for {}|{}."
            .format(codon_index, len(codons), taxon, gene)
        )

    codon_alignment = "".join(aligned_codons)

    if len(codon_alignment) % 3 != 0:
        raise ValueError(
            "Back-translated alignment length is not divisible by 3 for {}|{}."
            .format(taxon, gene)
        )

    for start in range(0, len(codon_alignment), 3):
        aligned_codon = codon_alignment[start:start + 3]

        if "-" in aligned_codon and aligned_codon != "---":
            raise ValueError(
                "Non-codon gap detected for {}|{} at alignment nt {}-{}."
                .format(taxon, gene, start + 1, start + 3)
            )

    return codon_alignment


taxa = read_taxa(taxa_file)
overlap_rows = coding_overlap_report(cds_bed, gene_order)
duplicated_reference_nt = sum(row["overlap_nt"] for row in overlap_rows)

temp_parent = outdir.parent
temp_dir = Path(tempfile.mkdtemp(
    prefix=outdir.name + ".tmp.",
    dir=str(temp_parent),
))

try:
    aa_input_dir = temp_dir / "protein_inputs"
    aa_aligned_dir = temp_dir / "protein_alignments"
    codon_aligned_dir = temp_dir / "codon_alignments"
    log_dir = temp_dir / "mafft_logs"

    for directory in (
        aa_input_dir,
        aa_aligned_dir,
        codon_aligned_dir,
        log_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    codon_alignments = {}
    gene_lengths = {}
    gene_commands = {}
    qc_rows = []

    for gene in gene_order:
        nucleotide_path = input_dir / (gene + ".fasta")
        nucleotide_by_taxon = read_gene_matrix(nucleotide_path, taxa)

        protein_records = []

        for taxon in taxa:
            sequence = nucleotide_by_taxon[taxon]
            stops = definite_stop_codons(sequence)

            if stops:
                raise ValueError(
                    "{}|{} contains definite stop codon(s) at amino-acid "
                    "positions {}."
                    .format(taxon, gene, ",".join(map(str, stops)))
                )

            protein_records.append(
                SeqRecord(
                    Seq(protein_for_alignment(sequence)),
                    id=taxon,
                    description="",
                )
            )

        aa_input = aa_input_dir / (gene + ".aa.unaligned.fasta")
        aa_aligned = aa_aligned_dir / (gene + ".aa.mafft.fasta")
        mafft_log = log_dir / (gene + ".mafft.log")

        written = SeqIO.write(protein_records, str(aa_input), "fasta")

        if written != len(taxa):
            raise ValueError(
                "Wrote {} protein records for {}; expected {}."
                .format(written, gene, len(taxa))
            )

        gene_commands[gene] = run_mafft(
            aa_input,
            aa_aligned,
            mafft_log,
        )

        aligned_records = list(SeqIO.parse(str(aa_aligned), "fasta"))
        aligned_by_taxon = {record.id: str(record.seq).upper()
                            for record in aligned_records}

        if len(aligned_by_taxon) != len(aligned_records):
            raise ValueError(
                "Duplicate identifiers in MAFFT output for {}."
                .format(gene)
            )

        if set(aligned_by_taxon) != set(taxa):
            missing = sorted(set(taxa) - set(aligned_by_taxon))
            extra = sorted(set(aligned_by_taxon) - set(taxa))

            raise ValueError(
                "MAFFT taxon mismatch for {}. Missing: {}; extra: {}."
                .format(
                    gene,
                    ",".join(missing) or "<none>",
                    ",".join(extra) or "<none>",
                )
            )

        aa_lengths = {len(aligned_by_taxon[taxon]) for taxon in taxa}

        if len(aa_lengths) != 1:
            raise ValueError(
                "Unequal protein alignment lengths for {}."
                .format(gene)
            )

        aa_length = aa_lengths.pop()

        for column in range(aa_length):
            if all(aligned_by_taxon[taxon][column] == "-" for taxon in taxa):
                raise ValueError(
                    "All-gap protein column {} in {}."
                    .format(column + 1, gene)
                )

        codon_records = []
        gene_codon_by_taxon = {}

        for taxon in taxa:
            codon_sequence = backtranslate(
                aligned_by_taxon[taxon],
                nucleotide_by_taxon[taxon],
                taxon,
                gene,
            )

            gene_codon_by_taxon[taxon] = codon_sequence
            codon_records.append(
                SeqRecord(
                    Seq(codon_sequence),
                    id=taxon,
                    description="",
                )
            )

            qc_rows.append({
                "gene": gene,
                "taxon": taxon,
                "unaligned_nt": len(nucleotide_by_taxon[taxon]),
                "aligned_nt": len(codon_sequence),
                "N_count": codon_sequence.count("N"),
                "gap_count": codon_sequence.count("-"),
                "percent_N": (
                    100.0 * codon_sequence.count("N") / len(codon_sequence)
                ),
                "percent_gap": (
                    100.0 * codon_sequence.count("-") / len(codon_sequence)
                ),
                "definite_internal_stop_count": 0,
            })

        codon_output = codon_aligned_dir / (gene + ".codon.mafft.fasta")
        written = SeqIO.write(codon_records, str(codon_output), "fasta")

        if written != len(taxa):
            raise ValueError(
                "Wrote {} codon records for {}; expected {}."
                .format(written, gene, len(taxa))
            )

        codon_lengths = {
            len(gene_codon_by_taxon[taxon])
            for taxon in taxa
        }

        if len(codon_lengths) != 1:
            raise ValueError(
                "Unequal codon alignment lengths for {}."
                .format(gene)
            )

        gene_length = codon_lengths.pop()

        if gene_length != aa_length * 3:
            raise ValueError(
                "{} codon alignment is {} nt but protein alignment is {} aa."
                .format(gene, gene_length, aa_length)
            )

        codon_alignments[gene] = gene_codon_by_taxon
        gene_lengths[gene] = gene_length

    # Concatenate completed codon alignments in a fixed biological order.
    concatenated_records = []

    for taxon in taxa:
        concatenated = "".join(
            codon_alignments[gene][taxon]
            for gene in gene_order
        )

        concatenated_records.append(
            SeqRecord(
                Seq(concatenated),
                id=taxon,
                description="",
            )
        )

    matrix_path = temp_dir / "Sturnira_13PCG_{}.codon_aligned.fasta".format(
        outdir.name
    )
    written = SeqIO.write(
        concatenated_records,
        str(matrix_path),
        "fasta",
    )

    if written != len(taxa):
        raise ValueError(
            "Wrote {} concatenated records; expected {}."
            .format(written, len(taxa))
        )

    matrix_records = list(SeqIO.parse(str(matrix_path), "fasta"))
    matrix_lengths = {len(record.seq) for record in matrix_records}

    if len(matrix_lengths) != 1:
        raise ValueError("Concatenated matrix has unequal sequence lengths.")

    matrix_length = matrix_lengths.pop()

    # Generate IQ-TREE gene-by-codon-position partitions.
    partition_path = temp_dir / "Sturnira_13PCG_{}.partitions.txt".format(
        outdir.name
    )
    coordinate_path = temp_dir / "Sturnira_13PCG_{}.coordinates.tsv".format(
        outdir.name
    )

    cursor = 1
    coordinate_rows = []

    with partition_path.open("w") as partition_handle:
        for gene in gene_order:
            start = cursor
            end = start + gene_lengths[gene] - 1

            partition_handle.write(
                "DNA, {0}_pos1 = {1}-{2}\\3\n".format(gene, start, end)
            )
            partition_handle.write(
                "DNA, {0}_pos2 = {1}-{2}\\3\n".format(gene, start + 1, end)
            )
            partition_handle.write(
                "DNA, {0}_pos3 = {1}-{2}\\3\n".format(gene, start + 2, end)
            )

            coordinate_rows.append({
                "gene": gene,
                "start_1based": start,
                "end_1based": end,
                "aligned_nt": gene_lengths[gene],
                "aligned_codons": gene_lengths[gene] // 3,
            })

            cursor = end + 1

    if cursor - 1 != matrix_length:
        raise ValueError(
            "Partition coordinates end at {}, matrix length is {}."
            .format(cursor - 1, matrix_length)
        )

    with coordinate_path.open("w", newline="") as handle:
        fieldnames = [
            "gene",
            "start_1based",
            "end_1based",
            "aligned_nt",
            "aligned_codons",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(coordinate_rows)

    qc_path = temp_dir / "Sturnira_13PCG_{}.alignment_qc.tsv".format(
        outdir.name
    )

    with qc_path.open("w", newline="") as handle:
        fieldnames = [
            "gene",
            "taxon",
            "unaligned_nt",
            "aligned_nt",
            "N_count",
            "gap_count",
            "percent_N",
            "percent_gap",
            "definite_internal_stop_count",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    overlap_path = temp_dir / "coding_overlap_report.tsv"

    with overlap_path.open("w", newline="") as handle:
        fieldnames = [
            "gene1",
            "gene2",
            "overlap_nt",
            "first_reference_position_1based",
            "last_reference_position_1based",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(overlap_rows)

    strategy_path = temp_dir / "alignment_strategy.tsv"

    with strategy_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "dataset",
            "taxa",
            "genes",
            "matrix_nt",
            "mafft_strategy",
            "translation_table",
            "full_gene_concatenation",
            "duplicated_reference_nt_from_overlaps",
        ])
        writer.writerow([
            outdir.name,
            len(taxa),
            len(gene_order),
            matrix_length,
            strategy,
            2,
            "yes",
            duplicated_reference_nt,
        ])

    command_path = temp_dir / "mafft_commands.tsv"

    with command_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "command"])

        for gene in gene_order:
            writer.writerow([gene, gene_commands[gene]])

    (temp_dir / "taxa.txt").write_text(
        "".join(taxon + "\n" for taxon in taxa)
    )

    # Replace the dataset output only after all genes and outputs validate.
    if outdir.exists():
        shutil.rmtree(str(outdir))

    shutil.move(str(temp_dir), str(outdir))

except Exception:
    shutil.rmtree(str(temp_dir), ignore_errors=True)
    raise

print("Dataset:             {}".format(outdir.name))
print("Taxa:                {}".format(len(taxa)))
print("Genes:               {}".format(len(gene_order)))
print("MAFFT strategy:      {}".format(strategy))
print("Concatenated length: {} nt".format(matrix_length))
print("Overlapping reference nt represented twice: {}".format(
    duplicated_reference_nt
))
print("Output directory:    {}".format(outdir))
PY

############################################
# FINAL CHECKS
############################################
MATRIX="$OUTDIR/Sturnira_13PCG_${LABEL}.codon_aligned.fasta"
PARTITIONS="$OUTDIR/Sturnira_13PCG_${LABEL}.partitions.txt"
COORDINATES="$OUTDIR/Sturnira_13PCG_${LABEL}.coordinates.tsv"
QC="$OUTDIR/Sturnira_13PCG_${LABEL}.alignment_qc.tsv"
STRATEGY_FILE="$OUTDIR/alignment_strategy.tsv"
OVERLAP_REPORT="$OUTDIR/coding_overlap_report.tsv"

for file in "$MATRIX" "$PARTITIONS" "$COORDINATES" "$QC" "$STRATEGY_FILE" "$OVERLAP_REPORT"; do
    [[ -s "$file" ]] || {
        echo "ERROR: expected output is missing or empty: $file" >&2
        exit 1
    }
done

N_EXPECTED=$(grep -cv '^[[:space:]]*$' "$TAXA_FILE" || true)
N_ALIGNED=$(grep -c '^>' "$MATRIX")

if [[ "$N_ALIGNED" -ne "$N_EXPECTED" ]]; then
    echo "ERROR: aligned matrix has $N_ALIGNED taxa; expected $N_EXPECTED." >&2
    exit 1
fi

N_PARTITIONS=$(grep -cv '^[[:space:]]*$' "$PARTITIONS" || true)

if [[ "$N_PARTITIONS" -ne 39 ]]; then
    echo "ERROR: partition file has $N_PARTITIONS partitions; expected 39." >&2
    exit 1
fi

echo
echo "$LABEL 13-PCG codon alignment completed successfully."
echo "Matrix:      $MATRIX"
echo "Partitions:  $PARTITIONS"
echo "QC report:   $QC"
echo "Overlap QC:  $OVERLAP_REPORT"
