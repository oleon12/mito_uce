#!/usr/bin/env bash
#SBATCH --job-name=mafft_complete_mt
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --array=0-3
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/mafft_complete_mt_%A_%a.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/mafft_complete_mt_%A_%a.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Align complete masked mitochondrial genomes for the optional M50, M40,
# M30, and M20 whole-mitogenome analyses.
#
# The sequences must already:
# - be in the same orientation;
# - be linearized at the same homologous circular-genome origin;
# - contain no pre-existing alignment gaps.
#
# This script adds the S. ludovici reference as a temporary alignment anchor,
# preserves an alignment containing that anchor for coordinate mapping, then
# removes the anchor from the final phylogenetic matrix.
#
# No automatic trimming is performed. The control region and other difficult
# regions must be inspected before downstream whole-mitogenome phylogenetic analysis.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

THRESHOLDS=(50 40 30 20)
THRESHOLD="${THRESHOLDS[$SLURM_ARRAY_TASK_ID]}"
LABEL="M${THRESHOLD}"

REF="references/S_ludovici_QCAZ_18312.fasta"

INPUT_ROOT="results/consensus_matrices"
INPUT="$INPUT_ROOT/$LABEL/Sturnira_complete_mt_${LABEL}.unaligned.fasta"
TAXA_FILE="$INPUT_ROOT/$LABEL/taxa.txt"

OUTROOT="results/alignments/consensus"
OUTDIR="$OUTROOT/$LABEL"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
REFERENCE_ALIGNMENT_ID="__STURNIRA_MT_REFERENCE__"

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
REF=$(realpath "$REF")
INPUT=$(realpath "$INPUT")
TAXA_FILE=$(realpath "$TAXA_FILE")
OUTROOT=$(realpath -m "$OUTROOT")
OUTDIR=$(realpath -m "$OUTDIR")

for file in "$REF" "$INPUT" "$TAXA_FILE"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$OUTROOT"

echo "Dataset:         $LABEL"
echo "Input:           $INPUT"
echo "Threads:         $THREADS"
echo "MAFFT strategy:  FFT-NS-i with 1,000 refinement iterations"

############################################
# ALIGN AND VALIDATE
############################################
python3 - \
    "$REF" \
    "$INPUT" \
    "$TAXA_FILE" \
    "$OUTDIR" \
    "$THREADS" \
    "$REFERENCE_ALIGNMENT_ID" <<'PY'
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


reference_path = Path(sys.argv[1])
input_path = Path(sys.argv[2])
taxa_file = Path(sys.argv[3])
outdir = Path(sys.argv[4])
threads = int(sys.argv[5])
reference_alignment_id = sys.argv[6]

allowed_nt = set("ACGTRYSWKMBDHVN")
unambiguous_nt = set("ACGT")


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )

    return records[0]


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


def validate_unaligned_sequence(sequence, taxon, source):
    if not sequence:
        raise ValueError(
            "Empty sequence for '{}' in {}."
            .format(taxon, source)
        )

    invalid = sorted(set(sequence) - allowed_nt)

    if invalid:
        raise ValueError(
            "Invalid or aligned symbols for '{}' in {}: {}."
            .format(taxon, source, ",".join(invalid))
        )

    if "-" in sequence:
        raise ValueError(
            "{} contains gaps; complete-genome inputs must be unaligned."
            .format(source)
        )


reference_record = read_single_fasta(reference_path)
reference_sequence = (
    str(reference_record.seq).upper().replace("U", "T")
)
validate_unaligned_sequence(
    reference_sequence,
    reference_record.id,
    reference_path,
)

taxa = read_taxa(taxa_file)
input_records = list(SeqIO.parse(str(input_path), "fasta"))

if not input_records:
    raise ValueError("No FASTA records found in {}.".format(input_path))

input_ids = [record.id for record in input_records]

if len(input_ids) != len(set(input_ids)):
    raise ValueError("Duplicate FASTA identifiers in {}.".format(input_path))

if input_ids != taxa:
    raise ValueError(
        "Taxon names/order in {} do not match {}."
        .format(input_path, taxa_file)
    )

if reference_alignment_id in set(taxa):
    raise ValueError(
        "Reserved reference alignment ID already occurs in taxa: {}."
        .format(reference_alignment_id)
    )

input_sequences = {}

for record in input_records:
    sequence = str(record.seq).upper().replace("U", "T")
    validate_unaligned_sequence(sequence, record.id, input_path)
    input_sequences[record.id] = sequence

temp_parent = outdir.parent
temp_dir = Path(tempfile.mkdtemp(
    prefix=outdir.name + ".tmp.",
    dir=str(temp_parent),
))

try:
    mafft_input = temp_dir / "complete_mt_with_reference.unaligned.fasta"
    mafft_output = temp_dir / "complete_mt_with_reference.mafft.fasta"
    mafft_log = temp_dir / "mafft.log"

    records_for_alignment = [
        SeqRecord(
            Seq(reference_sequence),
            id=reference_alignment_id,
            description="source={}".format(reference_record.id),
        )
    ] + [
        SeqRecord(
            Seq(input_sequences[taxon]),
            id=taxon,
            description="",
        )
        for taxon in taxa
    ]

    written = SeqIO.write(
        records_for_alignment,
        str(mafft_input),
        "fasta",
    )

    if written != len(taxa) + 1:
        raise ValueError(
            "Wrote {} MAFFT input records; expected {}."
            .format(written, len(taxa) + 1)
        )

    command = [
        "mafft",
        "--retree",
        "2",
        "--maxiterate",
        "1000",
        "--thread",
        str(threads),
        str(mafft_input),
    ]

    with mafft_output.open("w") as stdout, mafft_log.open("w") as stderr:
        completed = subprocess.run(
            command,
            stdout=stdout,
            stderr=stderr,
        )

    if completed.returncode != 0:
        raise RuntimeError(
            "MAFFT failed. See {}.".format(mafft_log)
        )

    aligned_records = list(SeqIO.parse(str(mafft_output), "fasta"))
    aligned_by_id = {
        record.id: str(record.seq).upper()
        for record in aligned_records
    }

    if len(aligned_by_id) != len(aligned_records):
        raise ValueError("Duplicate identifiers in MAFFT output.")

    expected_ids = [reference_alignment_id] + taxa

    if set(aligned_by_id) != set(expected_ids):
        missing = sorted(set(expected_ids) - set(aligned_by_id))
        extra = sorted(set(aligned_by_id) - set(expected_ids))

        raise ValueError(
            "MAFFT output ID mismatch. Missing: {}; extra: {}."
            .format(
                ",".join(missing) or "<none>",
                ",".join(extra) or "<none>",
            )
        )

    alignment_lengths = {len(sequence) for sequence in aligned_by_id.values()}

    if len(alignment_lengths) != 1:
        raise ValueError("MAFFT output sequences have unequal lengths.")

    alignment_length = alignment_lengths.pop()

    # MAFFT must preserve each ungapped input sequence exactly.
    if aligned_by_id[reference_alignment_id].replace("-", "") != reference_sequence:
        raise ValueError(
            "Ungapped aligned reference differs from input reference."
        )

    for taxon in taxa:
        if aligned_by_id[taxon].replace("-", "") != input_sequences[taxon]:
            raise ValueError(
                "Ungapped MAFFT output differs from input for {}."
                .format(taxon)
            )

    for column in range(alignment_length):
        if all(
            aligned_by_id[record_id][column] == "-"
            for record_id in expected_ids
        ):
            raise ValueError(
                "All-gap alignment column {}."
                .format(column + 1)
            )

    # Save the alignment including the temporary reference anchor.
    with_reference_path = (
        temp_dir
        / "Sturnira_complete_mt_{}.mafft.with_reference.fasta"
        .format(outdir.name)
    )

    with_reference_records = [
        SeqRecord(
            Seq(aligned_by_id[record_id]),
            id=record_id,
            description="",
        )
        for record_id in expected_ids
    ]

    SeqIO.write(
        with_reference_records,
        str(with_reference_path),
        "fasta",
    )

    # Final phylogenetic matrix excludes the temporary reference anchor.
    final_path = (
        temp_dir
        / "Sturnira_complete_mt_{}.mafft.fasta"
        .format(outdir.name)
    )

    final_records = [
        SeqRecord(
            Seq(aligned_by_id[taxon]),
            id=taxon,
            description="",
        )
        for taxon in taxa
    ]

    SeqIO.write(final_records, str(final_path), "fasta")

    # Map every reference genome position to an alignment column.
    coordinate_map_path = (
        temp_dir
        / "Sturnira_complete_mt_{}.reference_coordinate_map.tsv"
        .format(outdir.name)
    )

    reference_position = 0

    with coordinate_map_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "alignment_column_1based",
            "reference_position_1based",
            "reference_base",
        ])

        for column, base in enumerate(
            aligned_by_id[reference_alignment_id],
            start=1,
        ):
            if base == "-":
                writer.writerow([column, "NA", "-"])
            else:
                reference_position += 1
                writer.writerow([column, reference_position, base])

    if reference_position != len(reference_sequence):
        raise ValueError(
            "Reference coordinate map ended at {}, expected {}."
            .format(reference_position, len(reference_sequence))
        )

    # Per-taxon alignment QC.
    qc_path = (
        temp_dir
        / "Sturnira_complete_mt_{}.alignment_qc.tsv"
        .format(outdir.name)
    )
    qc_rows = []

    for taxon in taxa:
        sequence = aligned_by_id[taxon]
        leading_gaps = len(sequence) - len(sequence.lstrip("-"))
        trailing_gaps = len(sequence) - len(sequence.rstrip("-"))
        gap_count = sequence.count("-")
        n_count = sequence.count("N")
        other_ambiguity = sum(
            1 for base in sequence
            if base not in unambiguous_nt
            and base not in {"N", "-"}
        )

        flags = []

        if leading_gaps > 500:
            flags.append("long_leading_gap_check_circular_origin")

        if trailing_gaps > 500:
            flags.append("long_trailing_gap_check_circular_origin")

        if 100.0 * gap_count / alignment_length > 20.0:
            flags.append("high_gap_fraction")

        qc_rows.append({
            "taxon": taxon,
            "input_length_nt": len(input_sequences[taxon]),
            "alignment_length_nt": alignment_length,
            "gap_count": gap_count,
            "percent_gap": 100.0 * gap_count / alignment_length,
            "N_count": n_count,
            "percent_N_alignment": 100.0 * n_count / alignment_length,
            "other_ambiguity_count": other_ambiguity,
            "leading_gap_run": leading_gaps,
            "trailing_gap_run": trailing_gaps,
            "manual_review_flags": ";".join(flags) if flags else "none",
        })

    with qc_path.open("w", newline="") as handle:
        fieldnames = [
            "taxon",
            "input_length_nt",
            "alignment_length_nt",
            "gap_count",
            "percent_gap",
            "N_count",
            "percent_N_alignment",
            "other_ambiguity_count",
            "leading_gap_run",
            "trailing_gap_run",
            "manual_review_flags",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    # Per-column composition/missingness report for later inspection.
    column_stats_path = (
        temp_dir
        / "Sturnira_complete_mt_{}.column_stats.tsv"
        .format(outdir.name)
    )

    with column_stats_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "alignment_column_1based",
            "A",
            "C",
            "G",
            "T",
            "N",
            "gap",
            "other_ambiguity",
            "percent_gap",
            "percent_N",
        ])

        for column in range(alignment_length):
            characters = [
                aligned_by_id[taxon][column]
                for taxon in taxa
            ]

            counts = {
                base: characters.count(base)
                for base in ("A", "C", "G", "T", "N", "-")
            }
            other_ambiguity = sum(
                1 for base in characters
                if base not in {"A", "C", "G", "T", "N", "-"}
            )

            writer.writerow([
                column + 1,
                counts["A"],
                counts["C"],
                counts["G"],
                counts["T"],
                counts["N"],
                counts["-"],
                other_ambiguity,
                100.0 * counts["-"] / len(taxa),
                100.0 * counts["N"] / len(taxa),
            ])

    strategy_path = temp_dir / "alignment_strategy.tsv"

    with strategy_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "dataset",
            "taxa",
            "input_reference_length_nt",
            "alignment_length_nt",
            "mafft_strategy",
            "mafft_command",
            "automatic_trimming",
        ])
        writer.writerow([
            outdir.name,
            len(taxa),
            len(reference_sequence),
            alignment_length,
            "FFT-NS-i-1000",
            " ".join(command),
            "no",
        ])

    (temp_dir / "taxa.txt").write_text(
        "".join(taxon + "\n" for taxon in taxa)
    )

    if outdir.exists():
        shutil.rmtree(str(outdir))

    shutil.move(str(temp_dir), str(outdir))

except Exception:
    shutil.rmtree(str(temp_dir), ignore_errors=True)
    raise

flagged = [
    row["taxon"]
    for row in qc_rows
    if row["manual_review_flags"] != "none"
]

print("Dataset:              {}".format(outdir.name))
print("Taxa:                 {}".format(len(taxa)))
print("Reference length:     {} nt".format(len(reference_sequence)))
print("Alignment length:     {} nt".format(alignment_length))
print("Taxa flagged:         {}".format(len(flagged)))
print("Output directory:     {}".format(outdir))

if flagged:
    print(
        "WARNING: inspect manual_review_flags for: {}"
        .format(", ".join(flagged)),
        file=sys.stderr,
    )
PY

############################################
# FINAL CHECKS
############################################
MATRIX="$OUTDIR/Sturnira_complete_mt_${LABEL}.mafft.fasta"
WITH_REFERENCE="$OUTDIR/Sturnira_complete_mt_${LABEL}.mafft.with_reference.fasta"
COORDINATE_MAP="$OUTDIR/Sturnira_complete_mt_${LABEL}.reference_coordinate_map.tsv"
QC="$OUTDIR/Sturnira_complete_mt_${LABEL}.alignment_qc.tsv"
COLUMN_STATS="$OUTDIR/Sturnira_complete_mt_${LABEL}.column_stats.tsv"
STRATEGY_FILE="$OUTDIR/alignment_strategy.tsv"

for file in \
    "$MATRIX" \
    "$WITH_REFERENCE" \
    "$COORDINATE_MAP" \
    "$QC" \
    "$COLUMN_STATS" \
    "$STRATEGY_FILE"; do

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

echo
echo "$LABEL complete-mitogenome alignment completed."
echo "Matrix:           $MATRIX"
echo "QC report:        $QC"
echo "Coordinate map:   $COORDINATE_MAP"
echo
echo "Manual inspection is required before downstream whole-mitogenome phylogenetic analysis."
