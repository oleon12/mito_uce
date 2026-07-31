#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Independently validate every finalized 13-PCG FASTA for ingroups and
# outgroups. This test is read-only and writes:
#   results/tests/test_extract_cds.tsv

WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
OUTGROUP_LIST="results/outgroup_samples.txt"
CDS_METADATA="results/cds_coords.tsv"
INGROUP_SUMMARY="results/extract_masked_cds_summary.tsv"
OUTGROUP_SUMMARY="results/outgroup_extraction_summary.tsv"
TEST_DIR="results/tests"
TEST_REPORT="$TEST_DIR/test_extract_cds.tsv"

EXPECTED_GENES=(ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB)

cd "$WORKDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found in $CONDA_ENV." >&2
    exit 1
}

python3 - <<'PY_CHECK'
try:
    import Bio
except ImportError as exc:
    raise SystemExit("ERROR: Biopython is required.") from exc
PY_CHECK

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
OUTGROUP_LIST=$(realpath -m "$OUTGROUP_LIST")
CDS_METADATA=$(realpath "$CDS_METADATA")
INGROUP_SUMMARY=$(realpath "$INGROUP_SUMMARY")
OUTGROUP_SUMMARY=$(realpath -m "$OUTGROUP_SUMMARY")
TEST_DIR=$(realpath -m "$TEST_DIR")
TEST_REPORT=$(realpath -m "$TEST_REPORT")

for file in "$SAMPLE_LIST" "$CDS_METADATA" "$INGROUP_SUMMARY"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$TEST_DIR"
GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

python3 - \
    "$SAMPLE_LIST" \
    "$OUTGROUP_LIST" \
    "$CDS_METADATA" \
    "$INGROUP_SUMMARY" \
    "$OUTGROUP_SUMMARY" \
    "$TEST_REPORT" \
    "$GENES_TEXT" <<'PY'
import csv
import math
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


ingroup_list_path = Path(sys.argv[1])
outgroup_list_path = Path(sys.argv[2])
metadata_path = Path(sys.argv[3])
ingroup_summary_path = Path(sys.argv[4])
outgroup_summary_path = Path(sys.argv[5])
report_path = Path(sys.argv[6])
expected_genes = [x for x in sys.argv[7].split(",") if x]
results_dir = report_path.parent.parent


def read_list(path, three_fields=False):
    values = []
    seen = set()
    if not path.is_file():
        return values
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if three_fields and len(fields) != 3:
                raise ValueError(
                    "{} line {} must contain SAMPLE R1 R2."
                    .format(path, line_number)
                )
            value = fields[0]
            if value in seen:
                raise ValueError(
                    "Duplicate identifier '{}' in {}."
                    .format(value, path)
                )
            seen.add(value)
            values.append(value)
    return values


def read_summary(path, required, optional=False):
    if optional and not path.is_file():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "{} is missing columns: {}."
                .format(path, ", ".join(sorted(missing)))
            )
        rows = {}
        for row in reader:
            sample = row["sample"].strip()
            if sample in rows:
                raise ValueError(
                    "Duplicate summary row for '{}' in {}."
                    .format(sample, path)
                )
            rows[sample] = row
    return rows


def definite_stop_positions(sequence):
    positions = []
    for start in range(0, len(sequence), 3):
        codon = sequence[start:start + 3]
        if set(codon).issubset(set("ACGT")):
            if str(Seq(codon).translate(table=2)) == "*":
                positions.append(start // 3 + 1)
    return positions


expected_lengths = {}
with metadata_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"gene", "coding_nt_retained", "translation_validated"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            "{} is missing columns: {}."
            .format(metadata_path, ", ".join(sorted(missing)))
        )
    for row in reader:
        gene = row["gene"].strip().upper()
        if gene in expected_lengths:
            raise ValueError("Duplicate CDS metadata for {}.".format(gene))
        if row["translation_validated"].strip().lower() != "yes":
            raise ValueError("Reference CDS {} is not validated.".format(gene))
        expected_lengths[gene] = int(row["coding_nt_retained"])

if list(expected_lengths) != expected_genes:
    raise ValueError(
        "CDS metadata gene order/content differs from expected 13 PCGs."
    )

expected_total = sum(expected_lengths.values())
ingroups = read_list(ingroup_list_path, three_fields=True)
outgroups = read_list(outgroup_list_path)

if set(ingroups) & set(outgroups):
    raise ValueError("An identifier occurs as both ingroup and outgroup.")

ingroup_summary = read_summary(
    ingroup_summary_path,
    ["sample", "status", "total_coding_nt", "total_N", "percent_N_13PCG"],
)
outgroup_summary = read_summary(
    outgroup_summary_path,
    ["sample", "status", "total_13PCG_nt", "total_13PCG_N"],
    optional=True,
)

if set(ingroup_summary) != set(ingroups):
    raise ValueError("Ingroup summary and sample list differ.")
if outgroups and set(outgroup_summary) != set(outgroups):
    raise ValueError("Outgroup summary and outgroup list differ.")

rows = []
failures = []

for sample_type, samples in (("ingroup", ingroups), ("outgroup", outgroups)):
    for sample in samples:
        fasta = (
            results_dir / sample / "genes_masked"
            / (sample + "_13PCG_masked.fasta")
        )
        status = "failed"
        message = ""
        total_nt = ""
        total_n = ""
        percent_n = ""

        try:
            summary = (
                ingroup_summary[sample]
                if sample_type == "ingroup"
                else outgroup_summary[sample]
            )
            if summary["status"].strip().lower() != "passed":
                raise ValueError(
                    "Authoritative {} status is '{}': {}"
                    .format(
                        sample_type,
                        summary["status"],
                        summary.get("message", ""),
                    )
                )
            if not fasta.is_file():
                raise FileNotFoundError("Missing FASTA: {}.".format(fasta))

            records = list(SeqIO.parse(str(fasta), "fasta"))
            expected_ids = [
                "{}|{}".format(sample, gene) for gene in expected_genes
            ]
            ids = [record.id for record in records]

            if ids != expected_ids:
                raise ValueError(
                    "FASTA identifiers/order differ from the 13-gene contract."
                )

            sequences = {}
            for gene, record in zip(expected_genes, records):
                sequence = str(record.seq).upper().replace("U", "T")
                invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN"))
                if invalid:
                    raise ValueError(
                        "{} has invalid symbols: {}."
                        .format(gene, ",".join(invalid))
                    )
                if "-" in sequence:
                    raise ValueError("{} contains alignment gaps.".format(gene))
                if len(sequence) % 3 != 0:
                    raise ValueError(
                        "{} length {} is not divisible by 3."
                        .format(gene, len(sequence))
                    )
                stops = definite_stop_positions(sequence)
                if stops:
                    raise ValueError(
                        "{} has definite internal stop(s) at {}."
                        .format(gene, ",".join(map(str, stops)))
                    )
                if sample_type == "ingroup" and len(sequence) != expected_lengths[gene]:
                    raise ValueError(
                        "{} ingroup length {} differs from reference coding "
                        "length {}."
                        .format(gene, len(sequence), expected_lengths[gene])
                    )
                sequences[gene] = sequence

            total_nt = sum(len(sequence) for sequence in sequences.values())
            total_n = sum(sequence.count("N") for sequence in sequences.values())
            percent_n = 100.0 * total_n / total_nt

            if sample_type == "ingroup":
                if total_nt != expected_total:
                    raise ValueError(
                        "Ingroup total {} differs from expected {}."
                        .format(total_nt, expected_total)
                    )
                if int(summary["total_coding_nt"]) != total_nt:
                    raise ValueError("Ingroup summary coding length mismatch.")
                if int(summary["total_N"]) != total_n:
                    raise ValueError("Ingroup summary N count mismatch.")
                reported = float(summary["percent_N_13PCG"])
                if abs(reported - percent_n) > 1e-6:
                    raise ValueError("Ingroup summary percent N mismatch.")
            else:
                if int(summary["total_13PCG_nt"]) != total_nt:
                    raise ValueError("Outgroup summary coding length mismatch.")
                if int(summary["total_13PCG_N"]) != total_n:
                    raise ValueError("Outgroup summary N count mismatch.")

            status = "passed"
            message = "13-PCG structure, frame, stops, and summary agree"

        except Exception as exc:
            message = str(exc).replace("\t", " ").replace("\n", " ")
            failures.append(sample)

        rows.append({
            "sample": sample,
            "sample_type": sample_type,
            "status": status,
            "message": message,
            "fasta": str(fasta),
            "genes": 13 if status == "passed" else "",
            "total_nt": total_nt,
            "total_N": total_n,
            "percent_N": percent_n,
        })

with report_path.open("w", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("13-PCG FASTAs tested: {}".format(len(rows)))
print("Passed: {}".format(len(rows) - len(failures)))
print("Failed: {}".format(len(failures)))
print("Report: {}".format(report_path))

if failures:
    raise SystemExit(
        "ERROR: CDS test failures: {}."
        .format(", ".join(failures))
    )
PY

echo "test_extract_cds.sh passed."
