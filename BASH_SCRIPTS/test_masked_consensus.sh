#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Independently validate finalized ingroup masked consensuses and compare them
# with results/masked_consensus_summary.tsv. This script is read-only and writes:
#   results/tests/test_masked_consensus.tsv

WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
AUTHORITATIVE_SUMMARY="results/masked_consensus_summary.tsv"
TEST_DIR="results/tests"
TEST_REPORT="$TEST_DIR/test_masked_consensus.tsv"

cd "$WORKDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 samtools; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: $command_name was not found in $CONDA_ENV." >&2
        exit 1
    }
done

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
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
report_path = Path(sys.argv[4])
results_dir = report_path.parent.parent


def read_single(path):
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )
    return records[0]


def read_samples(path):
    values = []
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
            values.append(sample)
    return values


reference = read_single(reference_path)
reference_length = len(reference.seq)
samples = read_samples(sample_list_path)

with summary_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "sample", "status", "length_nt", "N_count", "percent_N",
        "combined_masked_positions", "consensus_file",
    }
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            "Masked-consensus summary is missing columns: {}."
            .format(", ".join(sorted(missing)))
        )
    summary = {row["sample"]: row for row in reader}

if set(summary) != set(samples):
    raise ValueError("Masked-consensus summary and sample list differ.")

rows = []
failures = []

for sample in samples:
    fasta = (
        results_dir / sample / "consensus_masked"
        / (sample + "_mt_masked.fasta")
    )
    status = "failed"
    message = ""
    length = n_count = percent_n = ""

    try:
        source = summary[sample]
        if source["status"].strip().lower() != "passed":
            raise ValueError(
                "Authoritative status is '{}': {}"
                .format(source["status"], source.get("message", ""))
            )
        if Path(source["consensus_file"]).resolve() != fasta.resolve():
            raise ValueError("Summary consensus path differs from expected path.")
        if not fasta.is_file():
            raise FileNotFoundError("Missing FASTA: {}.".format(fasta))
        fai = Path(str(fasta) + ".fai")
        if not fai.is_file():
            raise FileNotFoundError("Missing FASTA index: {}.".format(fai))

        record = read_single(fasta)
        sequence = str(record.seq).upper().replace("U", "T")
        if record.id != sample:
            raise ValueError(
                "FASTA ID '{}' does not match sample '{}'."
                .format(record.id, sample)
            )
        invalid = sorted(set(sequence) - set("ACGTN"))
        if invalid:
            raise ValueError(
                "Consensus contains invalid symbols: {}."
                .format(",".join(invalid))
            )
        if "-" in sequence:
            raise ValueError("Consensus contains alignment gaps.")
        if len(sequence) != reference_length:
            raise ValueError(
                "Consensus length {} differs from reference {}."
                .format(len(sequence), reference_length)
            )

        with fai.open() as handle:
            fai_rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
        if len(fai_rows) != 1:
            raise ValueError("FAI must contain exactly one record.")
        if fai_rows[0][0] != sample or int(fai_rows[0][1]) != reference_length:
            raise ValueError("FAI identifier or length does not match FASTA.")

        fetched = subprocess.run(
            ["samtools", "faidx", str(fasta), sample],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if fetched.returncode != 0:
            raise ValueError(
                "samtools faidx retrieval failed: {}"
                .format(fetched.stderr.strip())
            )
        fetched_sequence = "".join(
            line.strip()
            for line in fetched.stdout.splitlines()
            if line and not line.startswith(">")
        ).upper()
        if fetched_sequence != sequence:
            raise ValueError("FAI-retrieved sequence differs from FASTA.")

        length = len(sequence)
        n_count = sequence.count("N")
        percent_n = 100.0 * n_count / length

        if int(source["length_nt"]) != length:
            raise ValueError("Summary length does not match FASTA.")
        if int(source["N_count"]) != n_count:
            raise ValueError("Summary N_count does not match FASTA.")
        if abs(float(source["percent_N"]) - percent_n) > 1e-6:
            raise ValueError("Summary percent_N does not match FASTA.")
        if int(source["combined_masked_positions"]) != n_count:
            raise ValueError(
                "combined_masked_positions does not equal final N count."
            )

        status = "passed"
        message = "fixed-coordinate FASTA, FAI, and summary agree"

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failures.append(sample)

    rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "fasta": str(fasta),
        "length_nt": length,
        "N_count": n_count,
        "percent_N": percent_n,
        "reference_length_nt": reference_length,
    })

with report_path.open("w", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Masked consensuses tested: {}".format(len(rows)))
print("Passed: {}".format(len(rows) - len(failures)))
print("Failed: {}".format(len(failures)))
print("Report: {}".format(report_path))

if failures:
    raise SystemExit(
        "ERROR: masked-consensus test failures: {}."
        .format(", ".join(failures))
    )
PY

echo "test_masked_consensus.sh passed."
