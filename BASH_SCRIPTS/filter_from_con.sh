#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Create nested <=50%, <=40%, <=30%, and <=20% missing-data sample lists for the
# optional complete-mitochondrial-consensus analysis.
#
# This script consumes the validated global summary that will be produced by
# the revised make_masked_consensus.slurm. It filters INGROUP samples only.
# Outgroups are added later by concat_cons.sh.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"

# Required output contract for the future revised make_masked_consensus.slurm.
INPUT="results/masked_consensus_summary.tsv"

OUTDIR="results"
REPORT="$OUTDIR/filter_from_con_report.tsv"
COUNTS="$OUTDIR/filter_from_con_counts.tsv"

THRESHOLDS=(50 40 30 20)

############################################
# SETUP
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

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
INPUT=$(realpath "$INPUT")
OUTDIR=$(realpath -m "$OUTDIR")
REPORT=$(realpath -m "$REPORT")
COUNTS=$(realpath -m "$COUNTS")

for file in "$SAMPLE_LIST" "$REF" "$INPUT"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input is missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$OUTDIR"

THRESHOLDS_TEXT=$(IFS=,; echo "${THRESHOLDS[*]}")

############################################
# VALIDATE SUMMARY AND WRITE FILTER LISTS
############################################
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$INPUT" \
    "$OUTDIR" \
    "$REPORT" \
    "$COUNTS" \
    "$THRESHOLDS_TEXT" <<'PY'
import csv
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
input_path = Path(sys.argv[3])
outdir = Path(sys.argv[4])
report_path = Path(sys.argv[5])
counts_path = Path(sys.argv[6])
thresholds = [float(x) for x in sys.argv[7].split(",") if x]

if thresholds != [50.0, 40.0, 30.0, 20.0]:
    raise SystemExit(
        "ERROR: expected thresholds 50,40,30,20; observed {}."
        .format(",".join(map(str, thresholds)))
    )


def read_sample_ids(path):
    samples = []
    seen = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            sample = line.split()[0]

            if sample in seen:
                raise ValueError(
                    "Duplicate sample '{}' in {} at line {}."
                    .format(sample, path, line_number)
                )

            seen.add(sample)
            samples.append(sample)

    if not samples:
        raise ValueError("No sample identifiers found in {}.".format(path))

    return samples


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )

    return records[0]


def require_columns(reader, path, columns):
    observed = set(reader.fieldnames or [])
    missing = set(columns) - observed

    if missing:
        raise ValueError(
            "{} is missing required columns: {}."
            .format(path, ", ".join(sorted(missing)))
        )


def parse_float(value, label, sample):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(
            "{} for sample '{}' is not numeric: {!r}."
            .format(label, sample, value)
        )

    if not math.isfinite(number):
        raise ValueError(
            "{} for sample '{}' is not finite: {!r}."
            .format(label, sample, value)
        )

    return number


def parse_int(value, label, sample):
    number = parse_float(value, label, sample)

    if not number.is_integer():
        raise ValueError(
            "{} for sample '{}' is not an integer: {!r}."
            .format(label, sample, value)
        )

    return int(number)


sample_order = read_sample_ids(sample_list_path)
reference_record = read_single_fasta(reference_path)
reference_length = len(reference_record.seq)

if reference_length <= 0:
    raise ValueError("Reference sequence is empty.")

rows_by_sample = {}

with input_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    require_columns(
        reader,
        input_path,
        [
            "sample",
            "status",
            "length_nt",
            "N_count",
            "percent_N",
            "consensus_file",
        ],
    )

    for line_number, row in enumerate(reader, start=2):
        sample = row["sample"].strip()

        if not sample:
            raise ValueError(
                "Blank sample identifier in {} at line {}."
                .format(input_path, line_number)
            )

        if sample in rows_by_sample:
            raise ValueError(
                "Duplicate sample '{}' in {}."
                .format(sample, input_path)
            )

        rows_by_sample[sample] = row

listed = set(sample_order)
observed = set(rows_by_sample)

missing_samples = [sample for sample in sample_order if sample not in observed]
unexpected_samples = sorted(observed - listed)

if missing_samples:
    raise ValueError(
        "Samples from {} are absent from {}: {}."
        .format(
            sample_list_path,
            input_path,
            ", ".join(missing_samples),
        )
    )

if unexpected_samples:
    raise ValueError(
        "Unexpected samples occur in {} but not in {}: {}."
        .format(
            input_path,
            sample_list_path,
            ", ".join(unexpected_samples),
        )
    )

validated = []

for sample in sample_order:
    row = rows_by_sample[sample]
    status = row["status"].strip().lower()

    if status != "passed":
        message = row.get("message", "").strip()
        raise ValueError(
            "Sample '{}' has consensus status '{}'. Resolve the upstream "
            "failure before filtering. Message: {}"
            .format(sample, status or "<blank>", message or "<none>")
        )

    length_nt = parse_int(
        row["length_nt"],
        "length_nt",
        sample,
    )
    n_count = parse_int(
        row["N_count"],
        "N_count",
        sample,
    )
    reported_percent = parse_float(
        row["percent_N"],
        "percent_N",
        sample,
    )
    consensus_file = Path(row["consensus_file"])

    if length_nt != reference_length:
        raise ValueError(
            "Sample '{}' has consensus length {}, expected reference length {}."
            .format(sample, length_nt, reference_length)
        )

    if n_count < 0 or n_count > length_nt:
        raise ValueError(
            "Sample '{}' has invalid N_count={} for length_nt={}."
            .format(sample, n_count, length_nt)
        )

    recomputed_percent = 100.0 * n_count / length_nt

    if not (0.0 <= reported_percent <= 100.0):
        raise ValueError(
            "Sample '{}' has percent_N outside 0-100: {}."
            .format(sample, reported_percent)
        )

    if abs(reported_percent - recomputed_percent) > 1e-6:
        raise ValueError(
            "Sample '{}' has inconsistent missingness: reported {:.10f}%, "
            "recomputed {:.10f}%."
            .format(sample, reported_percent, recomputed_percent)
        )

    if not consensus_file.is_file():
        raise ValueError(
            "Consensus file recorded for sample '{}' does not exist: {}."
            .format(sample, consensus_file)
        )

    record = read_single_fasta(consensus_file)
    sequence = str(record.seq).upper().replace("U", "T")

    if record.id != sample:
        raise ValueError(
            "Consensus FASTA ID '{}' does not match sample '{}': {}."
            .format(record.id, sample, consensus_file)
        )

    if len(sequence) != length_nt:
        raise ValueError(
            "Consensus file length for '{}' is {}, summary reports {}."
            .format(sample, len(sequence), length_nt)
        )

    actual_n_count = sequence.count("N")

    if actual_n_count != n_count:
        raise ValueError(
            "Consensus file for '{}' contains {} Ns, summary reports {}."
            .format(sample, actual_n_count, n_count)
        )

    validated.append({
        "sample": sample,
        "status": "passed",
        "length_nt": length_nt,
        "N_count": n_count,
        "percent_N": recomputed_percent,
        "consensus_file": str(consensus_file),
    })

keep = {
    threshold: [
        row["sample"]
        for row in validated
        if row["percent_N"] <= threshold
    ]
    for threshold in thresholds
}

drop = {
    threshold: [
        row["sample"]
        for row in validated
        if row["percent_N"] > threshold
    ]
    for threshold in thresholds
}

if not set(keep[20.0]).issubset(set(keep[30.0])):
    raise ValueError(
        "The <=20% whole-consensus set is not a subset of the <=30% set."
    )

if not set(keep[30.0]).issubset(set(keep[40.0])):
    raise ValueError(
        "The <=30% whole-consensus set is not a subset of the <=40% set."
    )

if not set(keep[40.0]).issubset(set(keep[50.0])):
    raise ValueError(
        "The <=40% whole-consensus set is not a subset of the <=50% set."
    )

for threshold in thresholds:
    if set(keep[threshold]) & set(drop[threshold]):
        raise ValueError(
            "Keep/drop overlap detected at threshold {}%."
            .format(threshold)
        )

    if set(keep[threshold]) | set(drop[threshold]) != set(sample_order):
        raise ValueError(
            "Keep/drop lists do not account for every sample at threshold {}%."
            .format(threshold)
        )

temp_dir = Path(tempfile.mkdtemp(
    prefix=".filter_from_con.",
    dir=str(outdir),
))

try:
    for threshold in thresholds:
        threshold_int = int(threshold)

        (temp_dir / "keep_samples_cons_le{}.txt".format(threshold_int)).write_text(
            "".join(sample + "\n" for sample in keep[threshold])
        )

        (temp_dir / "drop_samples_cons_gt{}.txt".format(threshold_int)).write_text(
            "".join(sample + "\n" for sample in drop[threshold])
        )

    temp_report = temp_dir / report_path.name
    with temp_report.open("w", newline="") as handle:
        fieldnames = [
            "sample",
            "status",
            "length_nt",
            "N_count",
            "percent_N",
            "consensus_file",
            "keep_le50",
            "keep_le40",
            "keep_le30",
            "keep_le20",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()

        for row in validated:
            writer.writerow({
                **row,
                "keep_le50": "yes" if row["sample"] in keep[50.0] else "no",
                "keep_le40": "yes" if row["sample"] in keep[40.0] else "no",
                "keep_le30": "yes" if row["sample"] in keep[30.0] else "no",
                "keep_le20": "yes" if row["sample"] in keep[20.0] else "no",
            })

    temp_counts = temp_dir / counts_path.name
    with temp_counts.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "threshold_percent_N",
            "kept",
            "dropped",
            "total",
            "reference_length_nt",
        ])

        for threshold in thresholds:
            writer.writerow([
                int(threshold),
                len(keep[threshold]),
                len(drop[threshold]),
                len(sample_order),
                reference_length,
            ])

    output_names = [
        "keep_samples_cons_le50.txt",
        "keep_samples_cons_le40.txt",
        "keep_samples_cons_le30.txt",
        "keep_samples_cons_le20.txt",
        "drop_samples_cons_gt50.txt",
        "drop_samples_cons_gt40.txt",
        "drop_samples_cons_gt30.txt",
        "drop_samples_cons_gt20.txt",
        report_path.name,
        counts_path.name,
    ]

    for name in output_names:
        source = temp_dir / name
        destination = outdir / name

        if not source.exists():
            raise ValueError(
                "Expected temporary output was not created: {}."
                .format(source)
            )

        os.replace(str(source), str(destination))

finally:
    shutil.rmtree(str(temp_dir), ignore_errors=True)

print("Reference ID:              {}".format(reference_record.id))
print("Reference length:          {} bp".format(reference_length))
print("Validated samples:         {}".format(len(sample_order)))

for threshold in thresholds:
    print(
        "<={}% missing: {} kept; {} dropped"
        .format(
            int(threshold),
            len(keep[threshold]),
            len(drop[threshold]),
        )
    )

print("Detailed report:           {}".format(report_path))
print("Threshold counts:          {}".format(counts_path))
PY

############################################
# FINAL SHELL-LEVEL CHECKS
############################################
for threshold in "${THRESHOLDS[@]}"; do
    KEEP="$OUTDIR/keep_samples_cons_le${threshold}.txt"
    DROP="$OUTDIR/drop_samples_cons_gt${threshold}.txt"

    [[ -f "$KEEP" ]] || {
        echo "ERROR: expected keep list was not created: $KEEP" >&2
        exit 1
    }

    [[ -f "$DROP" ]] || {
        echo "ERROR: expected drop list was not created: $DROP" >&2
        exit 1
    }

    N_KEEP=$(grep -cv '^[[:space:]]*$' "$KEEP" || true)
    N_DROP=$(grep -cv '^[[:space:]]*$' "$DROP" || true)

    echo "Consensus <=${threshold}%: kept=$N_KEEP dropped=$N_DROP"
done

for file in "$REPORT" "$COUNTS"; do
    [[ -s "$file" ]] || {
        echo "ERROR: expected report is missing or empty: $file" >&2
        exit 1
    }
done

echo
echo "Whole-consensus filtering completed successfully."
echo "The M20, M30, M40, and M50 ingroup sets are validated and nested."
