#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Independently validate the final soft-filtered haploid mitochondrial VCFs
# and compare record counts with results/vcf_summary.tsv. This script is
# read-only and writes:
#   results/tests/test_vcf.tsv

WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
VCF_SUMMARY="results/vcf_summary.tsv"
TEST_DIR="results/tests"
TEST_REPORT="$TEST_DIR/test_vcf.tsv"

MIN_VARIANT_QUAL=20
MIN_DEPTH=3
MIN_ALT_DEPTH=3
MIN_ALT_FRACTION=0.80

cd "$WORKDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 bcftools; do
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
VCF_SUMMARY=$(realpath "$VCF_SUMMARY")
TEST_DIR=$(realpath -m "$TEST_DIR")
TEST_REPORT=$(realpath -m "$TEST_REPORT")

for file in "$SAMPLE_LIST" "$REF" "$VCF_SUMMARY"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$TEST_DIR"

python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$VCF_SUMMARY" \
    "$TEST_REPORT" \
    "$MIN_VARIANT_QUAL" \
    "$MIN_DEPTH" \
    "$MIN_ALT_DEPTH" \
    "$MIN_ALT_FRACTION" <<'PY'
import csv
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
report_path = Path(sys.argv[4])
min_qual = float(sys.argv[5])
min_depth = int(sys.argv[6])
min_alt_depth = int(sys.argv[7])
min_alt_fraction = float(sys.argv[8])
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


def parse_int(value):
    if value in ("", "."):
        return None
    return int(value)


def parse_ad(value):
    if value in ("", "."):
        return None
    pieces = value.split(",")
    if len(pieces) != 2 or any(piece in ("", ".") for piece in pieces):
        return None
    return [int(piece) for piece in pieces]


reference_records = list(SeqIO.parse(str(reference_path), "fasta"))
if len(reference_records) != 1:
    raise ValueError("Reference must contain exactly one FASTA record.")
reference_id = reference_records[0].id
reference_length = len(reference_records[0].seq)
samples = read_samples(sample_list_path)

with summary_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = set(reader.fieldnames or [])

    required = {
        "sample",
        "status",
        "vcf_file",
        "pass_haploid_alt_snp_records",
    }

    missing = required - fieldnames

    if missing:
        raise ValueError(
            "VCF summary is missing columns: {}."
            .format(", ".join(sorted(missing)))
        )

    if not (
        "total_variant_records" in fieldnames
        or "total_records" in fieldnames
    ):
        raise ValueError(
            "VCF summary must contain total_variant_records "
            "or total_records."
        )

    summary = {}

    for row in reader:
        sample = row["sample"].strip()

        if sample in summary:
            raise ValueError(
                "Duplicate VCF summary row for '{}'."
                .format(sample)
            )

        summary[sample] = row

if set(summary) != set(samples):
    raise ValueError("VCF summary and sample list differ.")

required_filter_ids = {
    "LowQual", "LowDepth", "MissingGT", "ReferenceGT",
    "NonHaploidGT", "NonSNP", "MissingAlleleDepth",
    "LowAltDepth", "LowAltFraction", "SnpNearIndel",
}

rows = []
failures = []

for sample in samples:
    vcf = results_dir / sample / "vcf" / (sample + ".vcf.gz")
    status = "failed"
    message = ""
    total_records = pass_records = filtered_records = ""

    try:
        source = summary[sample]
        if source["status"].strip().lower() != "passed":
            raise ValueError(
                "Authoritative VCF status is '{}': {}"
                .format(source["status"], source.get("message", ""))
            )
        if Path(source["vcf_file"]).resolve() != vcf.resolve():
            raise ValueError("Summary VCF path differs from expected path.")
        if not vcf.is_file():
            raise FileNotFoundError("Missing VCF: {}.".format(vcf))
        if not Path(str(vcf) + ".csi").is_file():
            raise FileNotFoundError("Missing CSI index for {}.".format(vcf))

        samples_in_vcf = [
            line.strip()
            for line in command_output(
                ["bcftools", "query", "-l", str(vcf)]
            ).splitlines()
            if line.strip()
        ]
        if samples_in_vcf != [sample]:
            raise ValueError(
                "VCF samples {} differ from expected [{}]."
                .format(samples_in_vcf, sample)
            )

        header = command_output(["bcftools", "view", "-h", str(vcf)])
        observed_filter_ids = set()
        for line in header.splitlines():
            if line.startswith("##FILTER=<ID="):
                observed_filter_ids.add(line.split("##FILTER=<ID=", 1)[1].split(",", 1)[0])
        missing_filters = required_filter_ids - observed_filter_ids
        if missing_filters:
            raise ValueError(
                "VCF header is missing FILTER definitions: {}."
                .format(", ".join(sorted(missing_filters)))
            )

        indexed_records = int(
            command_output(["bcftools", "index", "-n", str(vcf)]).strip()
        )
        query = command_output([
            "bcftools", "query",
            "-f", r"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%AD]\n",
            str(vcf),
        ])

        total_records = 0
        pass_records = 0
        filtered_records = 0

        for line_number, line in enumerate(query.splitlines(), start=1):
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    "Unexpected bcftools query fields at record {}: {}."
                    .format(line_number, len(fields))
                )
            chrom, pos_text, ref, alt, qual_text, filt, gt, dp_text, ad_text = fields
            pos = int(pos_text)
            if chrom != reference_id:
                raise ValueError("Unexpected contig '{}'.".format(chrom))
            if pos < 1 or pos > reference_length:
                raise ValueError("VCF position {} outside reference.".format(pos))
            if filt in ("", "."):
                raise ValueError("Record {} has an unset FILTER value.".format(pos))
            if gt not in (".", "0", "1"):
                raise ValueError(
                    "Record {} has non-haploid/unexpected GT '{}'."
                    .format(pos, gt)
                )

            total_records += 1
            if filt == "PASS":
                pass_records += 1
                alts = alt.split(",")
                if not (
                    len(ref) == 1
                    and len(alts) == 1
                    and len(alts[0]) == 1
                    and ref in "ACGT"
                    and alts[0] in "ACGT"
                ):
                    raise ValueError("PASS record {} is not a biallelic SNP.".format(pos))
                if gt != "1":
                    raise ValueError("PASS record {} does not have GT=1.".format(pos))
                if qual_text == "." or float(qual_text) < min_qual:
                    raise ValueError("PASS record {} fails QUAL threshold.".format(pos))
                dp = parse_int(dp_text)
                ad = parse_ad(ad_text)
                if dp is None or dp < min_depth:
                    raise ValueError("PASS record {} fails DP threshold.".format(pos))
                if ad is None:
                    raise ValueError("PASS record {} lacks valid AD.".format(pos))
                ref_depth, alt_depth = ad
                if alt_depth < min_alt_depth:
                    raise ValueError("PASS record {} fails ALT-depth threshold.".format(pos))
                if ref_depth + alt_depth <= 0:
                    raise ValueError("PASS record {} has zero AD total.".format(pos))
                alt_fraction = alt_depth / (ref_depth + alt_depth)
                if alt_fraction < min_alt_fraction:
                    raise ValueError("PASS record {} fails ALT-fraction threshold.".format(pos))
            else:
                filtered_records += 1

        if indexed_records != total_records:
            raise ValueError(
                "CSI record count {} differs from queried count {}."
                .format(indexed_records, total_records)
            )
        summary_total_text = (
            source.get("total_variant_records", "").strip()
        )

        if not summary_total_text:
            summary_total_text = (
                source.get("total_records", "").strip()
            )

        if not summary_total_text:
            raise ValueError(
                "Both total_variant_records and total_records "
                "are blank in vcf_summary.tsv."
            )

        try:
            summary_total_records = int(summary_total_text)
        except ValueError as exc:
            raise ValueError(
                "VCF summary total-record count is not an integer: {!r}."
                .format(summary_total_text)
            ) from exc

        if summary_total_records != total_records:
            raise ValueError(
                "VCF summary total-record count {} differs from "
                "the actual VCF count {}."
                .format(summary_total_records, total_records)
            )
        if int(source["pass_haploid_alt_snp_records"]) != pass_records:
            raise ValueError("VCF summary PASS-SNP count mismatch.")

        status = "passed"
        message = "VCF index, sample, haploid genotypes, PASS rules, and summary agree"

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failures.append(sample)

    rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "vcf": str(vcf),
        "total_variant_records": total_records,
        "pass_haploid_alt_snp_records": pass_records,
        "filtered_records": filtered_records,
    })

with report_path.open("w", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("VCFs tested: {}".format(len(rows)))
print("Passed: {}".format(len(rows) - len(failures)))
print("Failed: {}".format(len(failures)))
print("Report: {}".format(report_path))

if failures:
    raise SystemExit(
        "ERROR: VCF test failures: {}."
        .format(", ".join(failures))
    )
PY

echo "test_vcf.sh passed."
