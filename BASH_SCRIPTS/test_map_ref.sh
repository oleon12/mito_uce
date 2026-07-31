#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Independently validate finalized mapping BAMs and compare their core counts
# with results/mapping_summary.tsv. This script is read-only and writes:
#   results/tests/test_map_ref.tsv

WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
MAPPING_SUMMARY="results/mapping_summary.tsv"
TEST_DIR="results/tests"
TEST_REPORT="$TEST_DIR/test_map_ref.tsv"

cd "$WORKDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 samtools; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: $command_name was not found in $CONDA_ENV." >&2
        exit 1
    }
done

SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
MAPPING_SUMMARY=$(realpath "$MAPPING_SUMMARY")
TEST_DIR=$(realpath -m "$TEST_DIR")
TEST_REPORT=$(realpath -m "$TEST_REPORT")

for file in "$SAMPLE_LIST" "$REF" "$MAPPING_SUMMARY" "${REF}.fai"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$TEST_DIR"

python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$MAPPING_SUMMARY" \
    "$TEST_REPORT" <<'PY'
import csv
import subprocess
import sys
from pathlib import Path


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
report_path = Path(sys.argv[4])
results_dir = report_path.parent.parent


def command_output(command):
    completed = subprocess.run(
        list(map(str, command)),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed: {} | {}"
            .format(" ".join(map(str, command)), completed.stderr.strip())
        )
    return completed.stdout


def count_reads(bam, include=None, exclude=None, mapq=None):
    command = ["samtools", "view", "-c"]
    if include is not None:
        command.extend(["-f", str(include)])
    if exclude is not None:
        command.extend(["-F", str(exclude)])
    if mapq is not None:
        command.extend(["-q", str(mapq)])
    command.append(str(bam))
    return int(command_output(command).strip())


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
    return samples


samples = read_samples(sample_list_path)
with open(str(reference_path) + ".fai") as handle:
    fai_rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
if len(fai_rows) != 1:
    raise ValueError("Reference FAI must contain exactly one record.")
reference_id = fai_rows[0][0]
reference_length = int(fai_rows[0][1])

with summary_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "sample", "status", "primary_mapped_read_segments",
        "duplicate_read_segments", "nonduplicate_mapped_read_segments",
    }
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            "Mapping summary is missing columns: {}."
            .format(", ".join(sorted(missing)))
        )
    summary = {row["sample"]: row for row in reader}

if set(summary) != set(samples):
    raise ValueError("Mapping summary and sample list differ.")

rows = []
failures = []

for sample in samples:
    bam = results_dir / sample / "bam" / (sample + ".mapped.bam")
    status = "failed"
    message = ""
    mapped = duplicates = nonduplicates = mapq20_nonduplicates = ""

    try:
        source = summary[sample]
        if source["status"].strip().lower() != "passed":
            raise ValueError(
                "Authoritative mapping status is '{}': {}"
                .format(source["status"], source.get("message", ""))
            )
        if not bam.is_file():
            raise FileNotFoundError("Missing BAM: {}.".format(bam))
        if not Path(str(bam) + ".bai").is_file():
            raise FileNotFoundError("Missing BAM index for {}.".format(bam))

        quickcheck = subprocess.run(
            ["samtools", "quickcheck", "-v", str(bam)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if quickcheck.returncode != 0:
            raise ValueError(
                "samtools quickcheck failed: {} {}"
                .format(quickcheck.stdout.strip(), quickcheck.stderr.strip())
            )

        header = command_output(["samtools", "view", "-H", str(bam)])
        sm_values = set()
        rg_count = 0
        markdup_present = False
        for line in header.splitlines():
            fields = line.split("\t")
            if fields and fields[0] == "@RG":
                rg_count += 1
                for field in fields[1:]:
                    if field.startswith("SM:"):
                        sm_values.add(field[3:])
            if fields and fields[0] == "@PG" and "markdup" in line.lower():
                markdup_present = True
        if rg_count < 1 or sm_values != {sample}:
            raise ValueError(
                "Read-group sample identity is invalid: RG={}, SM={}."
                .format(rg_count, sorted(sm_values))
            )
        if not markdup_present:
            raise ValueError("BAM header does not document markdup.")

        idxstats = command_output(["samtools", "idxstats", str(bam)])
        reference_seen = False
        for line in idxstats.splitlines():
            if not line.strip():
                continue
            chrom, length_text, mapped_text, unmapped_text = line.split("\t")
            if chrom == "*":
                continue
            length = int(length_text)
            contig_mapped = int(mapped_text)
            if chrom == reference_id:
                reference_seen = True
                if length != reference_length:
                    raise ValueError("BAM reference length mismatch.")
            elif contig_mapped > 0:
                raise ValueError(
                    "Mapped reads occur on unexpected contig '{}'."
                    .format(chrom)
                )
        if not reference_seen:
            raise ValueError("Reference contig is absent from BAM.")

        unmapped = count_reads(bam, include=4)
        secondary = count_reads(bam, include=256)
        supplementary = count_reads(bam, include=2048)
        if unmapped or secondary or supplementary:
            raise ValueError(
                "Excluded alignment classes remain: unmapped={}, secondary={}, "
                "supplementary={}."
                .format(unmapped, secondary, supplementary)
            )

        mapped = count_reads(bam)
        duplicates = count_reads(bam, include=1024)
        nonduplicates = count_reads(bam, exclude=1024)
        mapq20_nonduplicates = count_reads(bam, exclude=1024, mapq=20)

        if mapped != duplicates + nonduplicates:
            raise ValueError("Duplicate/nonduplicate counts do not sum to BAM total.")
        if mapped != int(source["primary_mapped_read_segments"]):
            raise ValueError("Mapped-read count differs from mapping summary.")
        if duplicates != int(source["duplicate_read_segments"]):
            raise ValueError("Duplicate count differs from mapping summary.")
        if nonduplicates != int(source["nonduplicate_mapped_read_segments"]):
            raise ValueError("Nonduplicate count differs from mapping summary.")

        status = "passed"
        message = "BAM integrity, identity, flags, and summary counts agree"

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failures.append(sample)

    rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "bam": str(bam),
        "primary_mapped_read_segments": mapped,
        "duplicate_read_segments": duplicates,
        "nonduplicate_mapped_read_segments": nonduplicates,
        "mapq20_nonduplicate_read_segments": mapq20_nonduplicates,
    })

with report_path.open("w", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("BAMs tested: {}".format(len(rows)))
print("Passed: {}".format(len(rows) - len(failures)))
print("Failed: {}".format(len(failures)))
print("Report: {}".format(report_path))

if failures:
    raise SystemExit(
        "ERROR: mapping test failures: {}."
        .format(", ".join(failures))
    )
PY

echo "test_map_ref.sh passed."
