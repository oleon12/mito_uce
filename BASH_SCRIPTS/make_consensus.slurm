#!/usr/bin/env bash
#SBATCH --job-name=cons_mt_diagnostic
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/cons_mt_diagnostic_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/cons_mt_diagnostic_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Build an OPTIONAL, UNMASKED, SNP-only diagnostic consensus for each ingroup.
#
# DO NOT use these sequences for phylogenetic inference:
# positions with insufficient or absent coverage remain equal to the reference.
#
# The publication phylogenetic inputs are produced by:
#   make_masked_consensus.slurm
#
# This diagnostic script is retained to:
# - inspect which high-quality SNPs were applied;
# - compare masked and unmasked reconstructions;
# - detect variant/reference inconsistencies;
# - preserve a transparent intermediate product.
#
# To preserve reference coordinates, indels and other non-SNP variants are
# never applied here.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
RESULTS_DIR="results"

VCF_SUBDIR="vcf"
VCF_SUFFIX=".vcf.gz"

OUT_SUBDIR="consensus_diagnostic"
OUT_SUFFIX="_mt_unmasked_snp_only.fasta"

SUMMARY="$RESULTS_DIR/unmasked_consensus_summary.tsv"
FAILED_LIST="$RESULTS_DIR/unmasked_consensus_failed_samples.txt"

MIN_VARIANT_QUAL=20
REQUIRE_VCF_SAMPLE_MATCH="true"
OVERWRITE="true"
FAIL_JOB_ON_SAMPLE_ERROR="true"

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 bcftools samtools; do
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

echo "BCFtools version:"
bcftools --version | head -n 1
echo "SAMtools version:"
samtools --version | head -n 1

############################################
# RESOLVE AND CHECK PATHS
############################################
SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
RESULTS_DIR=$(realpath -m "$RESULTS_DIR")
SUMMARY=$(realpath -m "$SUMMARY")
FAILED_LIST=$(realpath -m "$FAILED_LIST")

for file in "$SAMPLE_LIST" "$REF"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input missing or empty: $file" >&2
        exit 1
    }
done

[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

mkdir -p "$RESULTS_DIR"

############################################
# BUILD DIAGNOSTIC CONSENSUSES
############################################
set +e
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$RESULTS_DIR" \
    "$VCF_SUBDIR" \
    "$VCF_SUFFIX" \
    "$OUT_SUBDIR" \
    "$OUT_SUFFIX" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$MIN_VARIANT_QUAL" \
    "$REQUIRE_VCF_SAMPLE_MATCH" \
    "$OVERWRITE" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" <<'PY'
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


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
results_dir = Path(sys.argv[3])
vcf_subdir = sys.argv[4]
vcf_suffix = sys.argv[5]
out_subdir = sys.argv[6]
out_suffix = sys.argv[7]
summary_path = Path(sys.argv[8])
failed_path = Path(sys.argv[9])
min_variant_qual = float(sys.argv[10])
require_vcf_sample_match = sys.argv[11].lower() == "true"
overwrite = sys.argv[12].lower() == "true"
fail_job_on_error = sys.argv[13].lower() == "true"

allowed_bases = set("ACGTRYSWKMBDHVN")


def run(command, log_handle, stdout_path=None):
    log_handle.write("$ {}\n".format(" ".join(map(str, command))))
    log_handle.flush()

    if stdout_path is None:
        completed = subprocess.run(
            list(map(str, command)),
            stdout=log_handle,
            stderr=log_handle,
        )
    else:
        with Path(stdout_path).open("w") as stdout:
            completed = subprocess.run(
                list(map(str, command)),
                stdout=stdout,
                stderr=log_handle,
            )

    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code {}: {}"
            .format(completed.returncode, " ".join(map(str, command)))
        )


def command_output(command):
    completed = subprocess.run(
        list(map(str, command)),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed: {}\n{}"
            .format(" ".join(map(str, command)), completed.stderr.strip())
        )

    return completed.stdout


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )

    return records[0]


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
        raise ValueError("No samples found in {}.".format(path))

    return samples


def count_vcf_records(path):
    output = command_output(["bcftools", "view", "-H", str(path)])
    return sum(1 for line in output.splitlines() if line.strip())


reference_record = read_single_fasta(reference_path)
reference_sequence = str(reference_record.seq).upper().replace("U", "T")
reference_length = len(reference_sequence)

if not reference_sequence:
    raise ValueError("Reference FASTA is empty.")

invalid_reference = sorted(set(reference_sequence) - allowed_bases)

if invalid_reference:
    raise ValueError(
        "Reference contains invalid symbols: {}."
        .format(",".join(invalid_reference))
    )

sample_ids = read_sample_ids(sample_list_path)
summary_rows = []
failed_samples = []

for sample in sample_ids:
    vcf_path = results_dir / sample / vcf_subdir / (sample + vcf_suffix)
    final_dir = results_dir / sample / out_subdir
    final_fasta = final_dir / (sample + out_suffix)

    status = "failed"
    message = ""
    vcf_sample = ""
    total_variant_records = ""
    hq_snp_records = ""
    final_length = ""
    final_n_count = ""

    if overwrite and final_dir.exists():
        shutil.rmtree(str(final_dir))

    temp_dir = None

    try:
        if not vcf_path.is_file():
            raise FileNotFoundError(
                "Missing VCF: {}".format(vcf_path)
            )

        vcf_samples = [
            line.strip()
            for line in command_output(
                ["bcftools", "query", "-l", str(vcf_path)]
            ).splitlines()
            if line.strip()
        ]

        if len(vcf_samples) != 1:
            raise ValueError(
                "{} must contain exactly one VCF sample; found {}."
                .format(vcf_path, len(vcf_samples))
            )

        vcf_sample = vcf_samples[0]

        if require_vcf_sample_match and vcf_sample != sample:
            raise ValueError(
                "VCF sample name '{}' does not match expected sample '{}'."
                .format(vcf_sample, sample)
            )

        sample_root = results_dir / sample
        sample_root.mkdir(parents=True, exist_ok=True)

        temp_dir = Path(tempfile.mkdtemp(
            prefix=".consensus_diagnostic.",
            dir=str(sample_root),
        ))

        normalized_vcf = temp_dir / (sample + ".normalized.vcf.gz")
        hq_snp_vcf = temp_dir / (sample + ".hq_haploid_snps.vcf.gz")
        raw_fasta = temp_dir / (sample + ".raw_reference_header.fasta")
        temp_final = temp_dir / final_fasta.name
        command_log = temp_dir / "commands.log"

        with command_log.open("w") as log:
            run([
                "bcftools", "norm",
                "-f", str(reference_path),
                "-c", "e",
                "-Oz",
                "-o", str(normalized_vcf),
                str(vcf_path),
            ], log)

            run([
                "bcftools", "index",
                "-f",
                str(normalized_vcf),
            ], log)

            # Retain only biallelic SNPs with:
            # - QUAL >= threshold;
            # - FILTER equal to PASS or unset (.);
            # - a haploid alternate genotype (GT="A").
            run([
                "bcftools", "view",
                "-s", vcf_sample,
                "-m2", "-M2",
                "-v", "snps",
                "-f", ".,PASS",
                "-i", 'QUAL>={} && GT="A"'.format(min_variant_qual),
                "-Oz",
                "-o", str(hq_snp_vcf),
                str(normalized_vcf),
            ], log)

            run([
                "bcftools", "index",
                "-f",
                str(hq_snp_vcf),
            ], log)

            run([
                "bcftools", "consensus",
                "-f", str(reference_path),
                "-s", vcf_sample,
                "-o", str(raw_fasta),
                str(hq_snp_vcf),
            ], log)

        total_variant_records = count_vcf_records(normalized_vcf)
        hq_snp_records = count_vcf_records(hq_snp_vcf)

        raw_record = read_single_fasta(raw_fasta)
        sequence = str(raw_record.seq).upper().replace("U", "T")

        if len(sequence) != reference_length:
            raise ValueError(
                "SNP-only consensus length {} differs from reference length {}."
                .format(len(sequence), reference_length)
            )

        invalid = sorted(set(sequence) - allowed_bases)

        if invalid:
            raise ValueError(
                "Consensus contains invalid symbols: {}."
                .format(",".join(invalid))
            )

        SeqIO.write(
            [SeqRecord(Seq(sequence), id=sample, description="")],
            str(temp_final),
            "fasta",
        )

        if overwrite:
            final_dir.mkdir(parents=True, exist_ok=True)

        # Preserve the normalized and filtered VCFs as auditable intermediates.
        final_dir.mkdir(parents=True, exist_ok=True)

        for source in (
            normalized_vcf,
            Path(str(normalized_vcf) + ".csi"),
            hq_snp_vcf,
            Path(str(hq_snp_vcf) + ".csi"),
            command_log,
            temp_final,
        ):
            if source.exists():
                os.replace(
                    str(source),
                    str(final_dir / source.name),
                )

        final_fasta = final_dir / temp_final.name
        final_length = len(sequence)
        final_n_count = sequence.count("N")
        status = "passed"
        message = (
            "DIAGNOSTIC_ONLY: low/no-depth positions remain reference; "
            "do not use for phylogenetic inference"
        )

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    finally:
        if temp_dir is not None:
            shutil.rmtree(str(temp_dir), ignore_errors=True)

    summary_rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "vcf_file": str(vcf_path),
        "vcf_sample": vcf_sample,
        "total_variant_records": total_variant_records,
        "hq_haploid_snp_records_applied": hq_snp_records,
        "length_nt": final_length,
        "N_count": final_n_count,
        "consensus_file": str(final_fasta),
        "phylogenetic_use": "no",
    })

summary_path.parent.mkdir(parents=True, exist_ok=True)

with summary_path.open("w", newline="") as handle:
    fieldnames = [
        "sample",
        "status",
        "message",
        "vcf_file",
        "vcf_sample",
        "total_variant_records",
        "hq_haploid_snp_records_applied",
        "length_nt",
        "N_count",
        "consensus_file",
        "phylogenetic_use",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

failed_path.write_text(
    "".join(sample + "\n" for sample in failed_samples)
)

passed_count = sum(row["status"] == "passed" for row in summary_rows)

print("Reference ID:          {}".format(reference_record.id))
print("Reference length:      {} bp".format(reference_length))
print("Samples listed:        {}".format(len(sample_ids)))
print("Samples passed:        {}".format(passed_count))
print("Samples failed:        {}".format(len(failed_samples)))
print("Summary:               {}".format(summary_path))
print("WARNING: these consensuses are diagnostic and not phylogenetic inputs.")

if failed_samples and fail_job_on_error:
    print(
        "ERROR: diagnostic consensus failures: {}."
        .format(", ".join(failed_samples)),
        file=sys.stderr,
    )
    sys.exit(2)
PY
PYTHON_STATUS=$?
set -e

[[ -s "$SUMMARY" ]] || {
    echo "ERROR: diagnostic summary missing or empty: $SUMMARY" >&2
    exit 1
}

[[ -f "$FAILED_LIST" ]] || {
    echo "ERROR: failed-sample list was not created: $FAILED_LIST" >&2
    exit 1
}

if [[ "$PYTHON_STATUS" -ne 0 ]]; then
    exit "$PYTHON_STATUS"
fi

echo
echo "Diagnostic SNP-only consensuses completed."
echo "Do not use them in MAFFT or IQ-TREE."
