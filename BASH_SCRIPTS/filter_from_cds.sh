#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Create nested <=50%, <=40%, <=30%, and <=20% missing-data sample lists for the
# primary 13-protein-coding-gene mitochondrial analysis.
#
# This script consumes the validated global summary produced by the revised
# extract_cds.slurm. It filters INGROUP samples only. Outgroups are added later
# by concat_cds.sh and are not included in these keep/drop lists.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"

SAMPLE_LIST="CONFS/sample_list.txt"
INPUT="results/extract_masked_cds_summary.tsv"
CDS_METADATA="results/cds_coords.tsv"

OUTDIR="results"
REPORT="$OUTDIR/filter_from_cds_report.tsv"
COUNTS="$OUTDIR/filter_from_cds_counts.tsv"

THRESHOLDS=(50 40 30 20)

############################################
# SETUP
############################################
cd "$WORKDIR"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found." >&2
    exit 1
}

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
INPUT=$(realpath "$INPUT")
CDS_METADATA=$(realpath "$CDS_METADATA")
OUTDIR=$(realpath -m "$OUTDIR")
REPORT=$(realpath -m "$REPORT")
COUNTS=$(realpath -m "$COUNTS")

for file in "$SAMPLE_LIST" "$INPUT" "$CDS_METADATA"; do
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
    "$INPUT" \
    "$CDS_METADATA" \
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


sample_list_path = Path(sys.argv[1])
input_path = Path(sys.argv[2])
metadata_path = Path(sys.argv[3])
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

# Derive the expected total retained coding length from the publication-grade
# CDS metadata rather than hard-coding 11,376 bp.
expected_genes = []
expected_total_nt = 0

with metadata_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    require_columns(
        reader,
        metadata_path,
        ["gene", "coding_nt_retained", "translation_validated"],
    )

    for row in reader:
        gene = row["gene"].strip().upper()

        if not gene:
            raise ValueError("Blank gene name in {}.".format(metadata_path))

        if gene in expected_genes:
            raise ValueError(
                "Duplicate gene '{}' in {}.".format(gene, metadata_path)
            )

        if row["translation_validated"].strip().lower() != "yes":
            raise ValueError(
                "Gene '{}' is not marked translation_validated=yes."
                .format(gene)
            )

        coding_nt = parse_int(
            row["coding_nt_retained"],
            "coding_nt_retained",
            gene,
        )

        if coding_nt <= 0 or coding_nt % 3 != 0:
            raise ValueError(
                "Gene '{}' has invalid coding length {}."
                .format(gene, coding_nt)
            )

        expected_genes.append(gene)
        expected_total_nt += coding_nt

if len(expected_genes) != 13:
    raise ValueError(
        "Expected 13 CDS metadata rows; found {}."
        .format(len(expected_genes))
    )

rows_by_sample = {}

with input_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    require_columns(
        reader,
        input_path,
        [
            "sample",
            "status",
            "total_coding_nt",
            "total_N",
            "percent_N_13PCG",
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
            "Sample '{}' has extraction status '{}'. Resolve the upstream "
            "failure before filtering. Message: {}"
            .format(sample, status or "<blank>", message or "<none>")
        )

    total_nt = parse_int(
        row["total_coding_nt"],
        "total_coding_nt",
        sample,
    )
    total_n = parse_int(
        row["total_N"],
        "total_N",
        sample,
    )
    reported_percent = parse_float(
        row["percent_N_13PCG"],
        "percent_N_13PCG",
        sample,
    )

    if total_nt != expected_total_nt:
        raise ValueError(
            "Sample '{}' has total_coding_nt={}, expected {} from {}."
            .format(sample, total_nt, expected_total_nt, metadata_path)
        )

    if total_n < 0 or total_n > total_nt:
        raise ValueError(
            "Sample '{}' has invalid total_N={} for total_coding_nt={}."
            .format(sample, total_n, total_nt)
        )

    recomputed_percent = 100.0 * total_n / total_nt

    if not (0.0 <= reported_percent <= 100.0):
        raise ValueError(
            "Sample '{}' has percent_N_13PCG outside 0-100: {}."
            .format(sample, reported_percent)
        )

    if abs(reported_percent - recomputed_percent) > 1e-6:
        raise ValueError(
            "Sample '{}' has inconsistent missingness: reported {:.10f}%, "
            "recomputed {:.10f}%."
            .format(sample, reported_percent, recomputed_percent)
        )

    validated.append({
        "sample": sample,
        "status": "passed",
        "total_coding_nt": total_nt,
        "total_N": total_n,
        "percent_N_13PCG": recomputed_percent,
    })

keep = {
    threshold: [
        row["sample"]
        for row in validated
        if row["percent_N_13PCG"] <= threshold
    ]
    for threshold in thresholds
}

drop = {
    threshold: [
        row["sample"]
        for row in validated
        if row["percent_N_13PCG"] > threshold
    ]
    for threshold in thresholds
}

# Explicitly validate nested sensitivity datasets.
if not set(keep[20.0]).issubset(set(keep[30.0])):
    raise ValueError("The <=20% CDS set is not a subset of the <=30% set.")

if not set(keep[30.0]).issubset(set(keep[40.0])):
    raise ValueError("The <=30% CDS set is not a subset of the <=40% set.")

if not set(keep[40.0]).issubset(set(keep[50.0])):
    raise ValueError("The <=40% CDS set is not a subset of the <=50% set.")

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

# Write all outputs to a temporary directory and replace them only after every
# validation has succeeded.
temp_dir = Path(tempfile.mkdtemp(
    prefix=".filter_from_cds.",
    dir=str(outdir),
))

try:
    for threshold in thresholds:
        threshold_int = int(threshold)

        (temp_dir / "keep_samples_cds_le{}.txt".format(threshold_int)).write_text(
            "".join(sample + "\n" for sample in keep[threshold])
        )

        (temp_dir / "drop_samples_cds_gt{}.txt".format(threshold_int)).write_text(
            "".join(sample + "\n" for sample in drop[threshold])
        )

    temp_report = temp_dir / report_path.name
    with temp_report.open("w", newline="") as handle:
        fieldnames = [
            "sample",
            "status",
            "total_coding_nt",
            "total_N",
            "percent_N_13PCG",
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
            "expected_total_coding_nt",
        ])

        for threshold in thresholds:
            writer.writerow([
                int(threshold),
                len(keep[threshold]),
                len(drop[threshold]),
                len(sample_order),
                expected_total_nt,
            ])

    output_names = [
        "keep_samples_cds_le50.txt",
        "keep_samples_cds_le40.txt",
        "keep_samples_cds_le30.txt",
        "keep_samples_cds_le20.txt",
        "drop_samples_cds_gt50.txt",
        "drop_samples_cds_gt40.txt",
        "drop_samples_cds_gt30.txt",
        "drop_samples_cds_gt20.txt",
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

print("Validated samples:          {}".format(len(sample_order)))
print("Expected 13-PCG length:     {} bp".format(expected_total_nt))

for threshold in thresholds:
    print(
        "<={}% missing: {} kept; {} dropped"
        .format(
            int(threshold),
            len(keep[threshold]),
            len(drop[threshold]),
        )
    )

print("Detailed report:            {}".format(report_path))
print("Threshold counts:           {}".format(counts_path))
PY

############################################
# FINAL SHELL-LEVEL CHECKS
############################################
for threshold in "${THRESHOLDS[@]}"; do
    KEEP="$OUTDIR/keep_samples_cds_le${threshold}.txt"
    DROP="$OUTDIR/drop_samples_cds_gt${threshold}.txt"

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

    echo "CDS <=${threshold}%: kept=$N_KEEP dropped=$N_DROP"
done

for file in "$REPORT" "$COUNTS"; do
    [[ -s "$file" ]] || {
        echo "ERROR: expected report is missing or empty: $file" >&2
        exit 1
    }
done

echo
echo "13-PCG filtering completed successfully."
echo "The M20, M30, M40, and M50 ingroup sets are validated and nested."
