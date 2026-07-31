#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Create matched-taxon intersections between the 13-PCG and complete-
# mitogenome missing-data filters at <=40%, <=30%, and <=20%.
#
# IMPORTANT:
# - These intersection lists are NOT used for the primary analyses.
# - Primary 13-PCG analyses use keep_samples_cds_le*.txt.
# - Primary complete-mitogenome analyses use keep_samples_cons_le*.txt.
# - Intersection lists are optional controls for comparing both data types
#   using exactly the same ingroup taxa.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"

SAMPLE_LIST="CONFS/sample_list.txt"
RESULTS_DIR="results"

THRESHOLDS=(40 30 20)

REPORT="$RESULTS_DIR/filter_intersection_report.tsv"
COUNTS="$RESULTS_DIR/filter_intersection_counts.tsv"

############################################
# SETUP
############################################
cd "$WORKDIR"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found." >&2
    exit 1
}

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
RESULTS_DIR=$(realpath -m "$RESULTS_DIR")
REPORT=$(realpath -m "$REPORT")
COUNTS=$(realpath -m "$COUNTS")

[[ -s "$SAMPLE_LIST" ]] || {
    echo "ERROR: sample list not found or empty: $SAMPLE_LIST" >&2
    exit 1
}

for threshold in "${THRESHOLDS[@]}"; do
    CDS_KEEP="$RESULTS_DIR/keep_samples_cds_le${threshold}.txt"
    CONS_KEEP="$RESULTS_DIR/keep_samples_cons_le${threshold}.txt"

    for file in "$CDS_KEEP" "$CONS_KEEP"; do
        [[ -f "$file" ]] || {
            echo "ERROR: required filter list not found: $file" >&2
            exit 1
        }
    done
done

mkdir -p "$RESULTS_DIR"

THRESHOLDS_TEXT=$(IFS=,; echo "${THRESHOLDS[*]}")

############################################
# VALIDATE AND BUILD INTERSECTIONS
############################################
python3 - \
    "$SAMPLE_LIST" \
    "$RESULTS_DIR" \
    "$REPORT" \
    "$COUNTS" \
    "$THRESHOLDS_TEXT" <<'PY'
import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path


sample_list_path = Path(sys.argv[1])
results_dir = Path(sys.argv[2])
report_path = Path(sys.argv[3])
counts_path = Path(sys.argv[4])
thresholds = [int(x) for x in sys.argv[5].split(",") if x]

if thresholds != [40, 30, 20]:
    raise SystemExit(
        "ERROR: expected thresholds 40,30,20; observed {}."
        .format(",".join(map(str, thresholds)))
    )


def read_ordered_unique(path, label):
    values = []
    seen = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            value = line.split()[0]

            if value in seen:
                raise ValueError(
                    "Duplicate sample '{}' in {} ({}) at line {}."
                    .format(value, path, label, line_number)
                )

            seen.add(value)
            values.append(value)

    return values


sample_order = read_ordered_unique(sample_list_path, "sample list")

if not sample_order:
    raise ValueError("No samples found in {}.".format(sample_list_path))

sample_set = set(sample_order)

cds = {}
cons = {}

for threshold in thresholds:
    cds_path = results_dir / "keep_samples_cds_le{}.txt".format(threshold)
    cons_path = results_dir / "keep_samples_cons_le{}.txt".format(threshold)

    cds_list = read_ordered_unique(cds_path, "CDS keep list")
    cons_list = read_ordered_unique(cons_path, "consensus keep list")

    unexpected_cds = sorted(set(cds_list) - sample_set)
    unexpected_cons = sorted(set(cons_list) - sample_set)

    if unexpected_cds:
        raise ValueError(
            "Unexpected samples in {}: {}."
            .format(cds_path, ", ".join(unexpected_cds))
        )

    if unexpected_cons:
        raise ValueError(
            "Unexpected samples in {}: {}."
            .format(cons_path, ", ".join(unexpected_cons))
        )

    cds[threshold] = set(cds_list)
    cons[threshold] = set(cons_list)

# Validate nestedness independently for both filtering branches.
if not cds[20].issubset(cds[30]) or not cds[30].issubset(cds[40]):
    raise ValueError("CDS keep lists are not nested: M20 ⊆ M30 ⊆ M40 failed.")

if not cons[20].issubset(cons[30]) or not cons[30].issubset(cons[40]):
    raise ValueError(
        "Whole-consensus keep lists are not nested: M20 ⊆ M30 ⊆ M40 failed."
    )

intersection = {
    threshold: cds[threshold] & cons[threshold]
    for threshold in thresholds
}

if not intersection[20].issubset(intersection[30]):
    raise ValueError(
        "Intersection M20 is not a subset of intersection M30."
    )

if not intersection[30].issubset(intersection[40]):
    raise ValueError(
        "Intersection M30 is not a subset of intersection M40."
    )

temp_dir = Path(tempfile.mkdtemp(
    prefix=".filter_intersection.",
    dir=str(results_dir),
))

try:
    report_rows = []

    for sample in sample_order:
        row = {"sample": sample}

        for threshold in thresholds:
            in_cds = sample in cds[threshold]
            in_cons = sample in cons[threshold]

            if in_cds and in_cons:
                category = "pass_both"
            elif in_cds:
                category = "pass_cds_only"
            elif in_cons:
                category = "pass_consensus_only"
            else:
                category = "fail_both"

            row["cds_le{}".format(threshold)] = "yes" if in_cds else "no"
            row["cons_le{}".format(threshold)] = "yes" if in_cons else "no"
            row["intersection_le{}".format(threshold)] = (
                "yes" if sample in intersection[threshold] else "no"
            )
            row["category_le{}".format(threshold)] = category

        report_rows.append(row)

    # Preserve the original sample-list order instead of sorting alphabetically.
    for threshold in thresholds:
        ordered_intersection = [
            sample for sample in sample_order
            if sample in intersection[threshold]
        ]

        output = (
            temp_dir
            / "keep_samples_intersection_le{}.txt".format(threshold)
        )
        output.write_text(
            "".join(sample + "\n" for sample in ordered_intersection)
        )

    report_temp = temp_dir / report_path.name
    report_fields = ["sample"]

    for threshold in thresholds:
        report_fields.extend([
            "cds_le{}".format(threshold),
            "cons_le{}".format(threshold),
            "intersection_le{}".format(threshold),
            "category_le{}".format(threshold),
        ])

    with report_temp.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=report_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    counts_temp = temp_dir / counts_path.name

    with counts_temp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "threshold_percent_N",
            "cds_kept",
            "consensus_kept",
            "intersection_kept",
            "cds_only",
            "consensus_only",
            "fail_both",
            "total_samples",
        ])

        for threshold in thresholds:
            cds_only = cds[threshold] - cons[threshold]
            cons_only = cons[threshold] - cds[threshold]
            fail_both = sample_set - (cds[threshold] | cons[threshold])

            writer.writerow([
                threshold,
                len(cds[threshold]),
                len(cons[threshold]),
                len(intersection[threshold]),
                len(cds_only),
                len(cons_only),
                len(fail_both),
                len(sample_order),
            ])

    output_names = [
        "keep_samples_intersection_le40.txt",
        "keep_samples_intersection_le30.txt",
        "keep_samples_intersection_le20.txt",
        report_path.name,
        counts_path.name,
    ]

    for name in output_names:
        source = temp_dir / name
        destination = results_dir / name

        if not source.exists():
            raise ValueError(
                "Expected temporary output was not created: {}."
                .format(source)
            )

        os.replace(str(source), str(destination))

finally:
    shutil.rmtree(str(temp_dir), ignore_errors=True)

for threshold in thresholds:
    print(
        "M{} matched intersection: {} samples"
        .format(threshold, len(intersection[threshold]))
    )

print("Detailed comparison: {}".format(report_path))
print("Counts summary:      {}".format(counts_path))
PY

############################################
# FINAL CHECKS
############################################
for threshold in "${THRESHOLDS[@]}"; do
    OUTPUT="$RESULTS_DIR/keep_samples_intersection_le${threshold}.txt"

    [[ -f "$OUTPUT" ]] || {
        echo "ERROR: expected intersection list missing: $OUTPUT" >&2
        exit 1
    }

    N=$(grep -cv '^[[:space:]]*$' "$OUTPUT" || true)
    echo "Intersection <=${threshold}%: $N samples"
done

for file in "$REPORT" "$COUNTS"; do
    [[ -s "$file" ]] || {
        echo "ERROR: expected report missing or empty: $file" >&2
        exit 1
    }
done

echo
echo "Matched-taxon intersection analysis completed."
echo "Do not use these lists for the primary CDS or complete-mitogenome trees."
echo "Use them only for optional same-taxon comparisons between data types."
