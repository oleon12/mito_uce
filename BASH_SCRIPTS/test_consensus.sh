#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Independently validate the OPTIONAL unmasked SNP-only diagnostic consensuses
# produced by make_consensus.slurm.
#
# These sequences are not phylogenetic inputs. This test never modifies them
# and writes only:
#   results/tests/test_consensus.tsv

WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
AUTHORITATIVE_SUMMARY="results/unmasked_consensus_summary.tsv"
TEST_DIR="results/tests"
TEST_REPORT="$TEST_DIR/test_consensus.tsv"

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
REF=$(realpath "$REF")
AUTHORITATIVE_SUMMARY=$(realpath "$AUTHORITATIVE_SUMMARY")
TEST_DIR=$(realpath -m "$TEST_DIR")
TEST_REPORT=$(realpath -m "$TEST_REPORT")

for file in "$SAMPLE_LIST" "$REF" "$AUTHORITATIVE_SUMMARY"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$TEST_DIR"

python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$AUTHORITATIVE_SUMMARY" \
    "$TEST_REPORT" <<'PY'
import csv
import math
import sys
from pathlib import Path

from Bio import SeqIO


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
report_path = Path(sys.argv[4])


def read_samples(path):
    samples = []
    seen = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(
                    "{} line {} must contain SAMPLE R1 R2."
                    .format(path, line_number)
                )
            sample = fields[0]
            if sample in seen:
                raise ValueError("Duplicate sample '{}'.".format(sample))
            seen.add(sample)
            samples.append(sample)
    if not samples:
        raise ValueError("No samples found in {}.".format(path))
    return samples


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )
    return records[0]


reference_record = read_single_fasta(reference_path)
reference_length = len(reference_record.seq)
samples = read_samples(sample_list_path)

with summary_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "sample", "status", "length_nt", "N_count",
        "consensus_file", "phylogenetic_use",
    }
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            "Authoritative summary is missing columns: {}."
            .format(", ".join(sorted(missing)))
        )
    summary = {}
    for row in reader:
        sample = row["sample"].strip()
        if sample in summary:
            raise ValueError("Duplicate summary row for '{}'.".format(sample))
        summary[sample] = row

if set(summary) != set(samples):
    raise ValueError(
        "Sample-list and unmasked-summary sample sets differ."
    )

rows = []
failures = []

for sample in samples:
    status = "failed"
    message = ""
    length = ""
    n_count = ""
    fasta_path = Path(summary[sample]["consensus_file"])

    try:
        source = summary[sample]
        if source["status"].strip().lower() != "passed":
            raise ValueError(
                "Authoritative status is '{}': {}"
                .format(source["status"], source.get("message", ""))
            )
        if source["phylogenetic_use"].strip().lower() != "no":
            raise ValueError(
                "Diagnostic consensus is not explicitly marked "
                "phylogenetic_use=no."
            )
        expected_path = (
            reference_path.parent.parent
            / "results" / sample / "consensus_diagnostic"
            / (sample + "_mt_unmasked_snp_only.fasta")
        )
        if fasta_path.resolve() != expected_path.resolve():
            raise ValueError(
                "Summary consensus path differs from expected path: {}."
                .format(expected_path)
            )
        if not fasta_path.is_file():
            raise FileNotFoundError("Missing FASTA: {}.".format(fasta_path))

        record = read_single_fasta(fasta_path)
        sequence = str(record.seq).upper().replace("U", "T")

        if record.id != sample:
            raise ValueError(
                "FASTA ID '{}' does not match sample '{}'."
                .format(record.id, sample)
            )
        if set(sequence) - set("ACGTN"):
            raise ValueError("FASTA contains invalid symbols.")
        if "-" in sequence:
            raise ValueError("Diagnostic consensus contains alignment gaps.")
        if len(sequence) != reference_length:
            raise ValueError(
                "Length {} differs from reference length {}."
                .format(len(sequence), reference_length)
            )

        length = len(sequence)
        n_count = sequence.count("N")

        if int(source["length_nt"]) != length:
            raise ValueError("Summary length does not match FASTA.")
        if int(source["N_count"]) != n_count:
            raise ValueError("Summary N_count does not match FASTA.")

        status = "passed"
        message = "diagnostic consensus structure and summary agree"

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failures.append(sample)

    rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "fasta": str(fasta_path),
        "length_nt": length,
        "N_count": n_count,
        "reference_length_nt": reference_length,
        "phylogenetic_use": "no",
    })

with report_path.open("w", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Diagnostic consensuses tested: {}".format(len(rows)))
print("Passed: {}".format(len(rows) - len(failures)))
print("Failed: {}".format(len(failures)))
print("Report: {}".format(report_path))

if failures:
    raise SystemExit(
        "ERROR: diagnostic consensus test failures: {}."
        .format(", ".join(failures))
    )
PY

echo "test_consensus.sh passed. Diagnostic sequences remain non-phylogenetic."
