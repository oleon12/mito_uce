#!/usr/bin/env bash
#SBATCH --job-name=vcf_mt_all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/vcf_mt_all_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/vcf_mt_all_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Call haploid mitochondrial SNPs and short indels from each publication-format
# mapped BAM and create a normalized, indexed, soft-filtered audit VCF.
#
# Final downstream VCF:
#   results/SAMPLE/vcf/SAMPLE.vcf.gz
#
# IMPORTANT:
# - The final VCF retains low-confidence and non-SNP records with FILTER labels.
# - PASS records are high-confidence haploid ALT SNP candidates.
# - make_masked_consensus.slurm applies only PASS haploid SNPs and masks every
#   other VCF record, preserving reference coordinates.
# - This is not a heteroplasmy caller. Allelic imbalance is used as a
#   conservative homoplasmic-consensus filter.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
RESULTS_DIR="results"

BAM_SUBDIR="bam"
BAM_SUFFIX=".mapped.bam"

VCF_SUBDIR="vcf"
FINAL_VCF_SUFFIX=".vcf.gz"

SUMMARY="$RESULTS_DIR/vcf_summary.tsv"
FAILED_LIST="$RESULTS_DIR/vcf_failed_samples.txt"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

# These evidence thresholds must remain synchronized with the publication
# masked-consensus protocol.
MIN_MAP_QUAL=20
MIN_BASE_QUAL=20
MIN_VARIANT_QUAL=20
MIN_DEPTH=3
MIN_ALT_DEPTH=3
MIN_ALT_FRACTION=0.80

# SNPs this many bases from a called indel/complex record are soft-filtered.
SNP_GAP_AROUND_NON_SNP=3

# Increase mpileup's default depth limit dynamically above the observed
# quality-filtered depth. The floor avoids the documented default cap of 250.
MPILEUP_DEPTH_FLOOR=1000
MPILEUP_DEPTH_MULTIPLIER=1.25
MPILEUP_DEPTH_PADDING=100

# Include anomalous/orphan read pairs in mpileup. This matches the evidence
# counted by samtools depth in the revised mapping and consensus scripts.
COUNT_ORPHANS="true"

# Full BAQ is conservative around local misalignments and makes the behavior
# explicit rather than relying on version-dependent heuristic application.
FULL_BAQ="true"

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

BCFTOOLS_VERSION="$(bcftools --version | awk 'NR==1 {print $2; exit}')"
SAMTOOLS_VERSION="$(samtools --version | awk 'NR==1 {print $2; exit}')"

echo "BCFtools version: $BCFTOOLS_VERSION"
echo "SAMtools version: $SAMTOOLS_VERSION"

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

[[ -s "${REF}.fai" ]] || {
    echo "ERROR: reference .fai missing. Run prep_ref.slurm first." >&2
    exit 1
}

mkdir -p "$RESULTS_DIR"

############################################
# CALL, NORMALIZE, SOFT-FILTER, AND VALIDATE
############################################
set +e
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$RESULTS_DIR" \
    "$BAM_SUBDIR" \
    "$BAM_SUFFIX" \
    "$VCF_SUBDIR" \
    "$FINAL_VCF_SUFFIX" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$THREADS" \
    "$MIN_MAP_QUAL" \
    "$MIN_BASE_QUAL" \
    "$MIN_VARIANT_QUAL" \
    "$MIN_DEPTH" \
    "$MIN_ALT_DEPTH" \
    "$MIN_ALT_FRACTION" \
    "$SNP_GAP_AROUND_NON_SNP" \
    "$MPILEUP_DEPTH_FLOOR" \
    "$MPILEUP_DEPTH_MULTIPLIER" \
    "$MPILEUP_DEPTH_PADDING" \
    "$COUNT_ORPHANS" \
    "$FULL_BAQ" \
    "$OVERWRITE" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" \
    "$BCFTOOLS_VERSION" \
    "$SAMTOOLS_VERSION" <<'PY'
import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from Bio import SeqIO


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
results_dir = Path(sys.argv[3])
bam_subdir = sys.argv[4]
bam_suffix = sys.argv[5]
vcf_subdir = sys.argv[6]
final_vcf_suffix = sys.argv[7]
summary_path = Path(sys.argv[8])
failed_path = Path(sys.argv[9])
threads = int(sys.argv[10])
min_map_qual = int(sys.argv[11])
min_base_qual = int(sys.argv[12])
min_variant_qual = float(sys.argv[13])
min_depth = int(sys.argv[14])
min_alt_depth = int(sys.argv[15])
min_alt_fraction = float(sys.argv[16])
snp_gap = int(sys.argv[17])
mpileup_depth_floor = int(sys.argv[18])
mpileup_depth_multiplier = float(sys.argv[19])
mpileup_depth_padding = int(sys.argv[20])
count_orphans = sys.argv[21].lower() == "true"
full_baq = sys.argv[22].lower() == "true"
overwrite = sys.argv[23].lower() == "true"
fail_job_on_error = sys.argv[24].lower() == "true"
bcftools_version = sys.argv[25]
samtools_version = sys.argv[26]

if threads < 1:
    raise SystemExit("ERROR: THREADS must be >=1.")

if min_depth < 1 or min_alt_depth < 1:
    raise SystemExit("ERROR: depth thresholds must be positive.")

if not 0.0 < min_alt_fraction <= 1.0:
    raise SystemExit("ERROR: MIN_ALT_FRACTION must lie in (0,1].")

if snp_gap < 0:
    raise SystemExit("ERROR: SNP_GAP_AROUND_NON_SNP must be >=0.")


FILTER_DEFINITIONS = {
    "LowQual": (
        "Variant QUAL is missing or below {:.3f}"
        .format(min_variant_qual)
    ),
    "LowDepth": (
        "Sample FORMAT/DP is missing or below {}"
        .format(min_depth)
    ),
    "MissingGT": "Sample genotype is missing",
    "ReferenceGT": "Variant record has a haploid reference genotype",
    "NonHaploidGT": "Genotype is not a single haploid allele",
    "NonSNP": "Record is not a biallelic single-nucleotide substitution",
    "MissingAlleleDepth": "FORMAT/AD is missing or malformed",
    "LowAltDepth": (
        "Alternate allele depth is below {}"
        .format(min_alt_depth)
    ),
    "LowAltFraction": (
        "Alternate allele fraction among AD-counted alleles is below {:.3f}"
        .format(min_alt_fraction)
    ),
    "SnpNearIndel": (
        "SNP lies within {} bp of an indel or complex variant"
        .format(snp_gap)
    ),
}


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


def run_command(command, log_handle, stdout_path=None):
    command = list(map(str, command))
    log_handle.write("$ {}\n".format(" ".join(command)))
    log_handle.flush()

    if stdout_path is None:
        completed = subprocess.run(
            command,
            stdout=log_handle,
            stderr=log_handle,
        )
    else:
        with Path(stdout_path).open("w") as stdout:
            completed = subprocess.run(
                command,
                stdout=stdout,
                stderr=log_handle,
            )

    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code {}: {}"
            .format(completed.returncode, " ".join(command))
        )


def run_pipeline(commands, log_handle):
    processes = []
    previous_stdout = None

    log_handle.write(
        "$ {}\n".format(
            " | ".join(" ".join(map(str, command)) for command in commands)
        )
    )
    log_handle.flush()

    try:
        for index, command in enumerate(commands):
            process = subprocess.Popen(
                list(map(str, command)),
                stdin=previous_stdout,
                stdout=(
                    subprocess.PIPE
                    if index < len(commands) - 1
                    else log_handle
                ),
                stderr=log_handle,
            )

            if previous_stdout is not None:
                previous_stdout.close()

            previous_stdout = process.stdout
            processes.append(process)

        return_codes = [process.wait() for process in processes]

    finally:
        if previous_stdout is not None:
            previous_stdout.close()

    failures = [
        (commands[index], return_code)
        for index, return_code in enumerate(return_codes)
        if return_code != 0
    ]

    if failures:
        details = "; ".join(
            "{} -> {}".format(" ".join(map(str, command)), return_code)
            for command, return_code in failures
        )
        raise RuntimeError("Pipeline failed: {}".format(details))


def read_sample_ids(path):
    samples = []
    seen = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()

            if len(fields) != 3:
                raise ValueError(
                    "{} line {} must contain exactly SAMPLE R1 R2."
                    .format(path, line_number)
                )

            sample = fields[0]

            if not re.fullmatch(r"[A-Za-z0-9._-]+", sample):
                raise ValueError(
                    "Unsafe sample identifier '{}'."
                    .format(sample)
                )

            if sample in seen:
                raise ValueError(
                    "Duplicate sample '{}' in {}."
                    .format(sample, path)
                )

            seen.add(sample)
            samples.append(sample)

    if not samples:
        raise ValueError("No samples found in {}.".format(path))

    return samples


def read_reference(path):
    records = list(SeqIO.parse(str(path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )

    record = records[0]
    sequence = str(record.seq).upper().replace("U", "T")

    if not sequence:
        raise ValueError("Reference sequence is empty.")

    if set(sequence) - set("ACGT"):
        raise ValueError("Reference contains non-ACGT symbols.")

    return record.id, len(sequence)


def parse_bam_header_sample(header_text):
    sample_names = set()
    read_group_count = 0
    markdup_present = False

    for line in header_text.splitlines():
        fields = line.split("\t")

        if not fields:
            continue

        if fields[0] == "@RG":
            read_group_count += 1

            for field in fields[1:]:
                if field.startswith("SM:"):
                    sample_names.add(field[3:])

        elif fields[0] == "@PG":
            lower = line.lower()

            if "markdup" in lower:
                markdup_present = True

    return read_group_count, sample_names, markdup_present


def observed_depths(bam, reference_id):
    output = command_output([
        "samtools", "depth",
        "-aa",
        "-s",
        "-q", str(min_base_qual),
        "-Q", str(min_map_qual),
        "-r", reference_id,
        str(bam),
    ])

    depths = []

    for line_number, line in enumerate(output.splitlines(), start=1):
        if not line.strip():
            continue

        fields = line.split("\t")

        if len(fields) != 3:
            raise ValueError(
                "Unexpected samtools depth output at line {}."
                .format(line_number)
            )

        chrom, position, depth = fields

        if chrom != reference_id:
            raise ValueError(
                "Unexpected depth contig '{}'."
                .format(chrom)
            )

        depths.append(int(depth))

    return depths


def parse_info(text):
    result = {}

    if text == "." or not text:
        return result

    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = True

    return result


def parse_format(format_text, sample_text):
    keys = format_text.split(":") if format_text else []
    values = sample_text.split(":") if sample_text else []

    if len(values) < len(keys):
        values.extend(["."] * (len(keys) - len(values)))

    return dict(zip(keys, values))


def parse_int(value):
    if value in (None, "", "."):
        return None

    return int(value)


def parse_int_list(value):
    if value in (None, "", "."):
        return None

    values = value.split(",")

    if any(item in ("", ".") for item in values):
        return None

    return [int(item) for item in values]


def variant_type(ref, alts):
    if (
        len(alts) == 1
        and len(ref) == 1
        and len(alts[0]) == 1
        and ref in "ACGT"
        and alts[0] in "ACGT"
    ):
        return "SNP"

    if any(len(alt) != len(ref) for alt in alts):
        return "INDEL"

    return "COMPLEX"


def circular_positions_around_interval(start_1based, ref_length, genome_length, flank):
    positions = set()
    end_1based = start_1based + ref_length - 1

    for position in range(start_1based - flank, end_1based + flank + 1):
        wrapped = ((position - 1) % genome_length) + 1
        positions.add(wrapped)

    return positions


def read_vcf_records(path):
    header_lines = []
    column_header = None
    records = []

    with gzip.open(str(path), "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if line.startswith("##"):
                header_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                column_header = line
                continue

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != 10:
                raise ValueError(
                    "{} line {} has {} columns; expected single-sample VCF "
                    "with 10 columns."
                    .format(path, line_number, len(fields))
                )

            chrom, pos, identifier, ref, alt, qual, filt, info, fmt, sample = fields
            alts = alt.split(",")

            records.append({
                "fields": fields,
                "chrom": chrom,
                "pos": int(pos),
                "id": identifier,
                "ref": ref,
                "alts": alts,
                "qual_text": qual,
                "filter_text": filt,
                "info": parse_info(info),
                "format": parse_format(fmt, sample),
                "type": variant_type(ref, alts),
            })

    if column_header is None:
        raise ValueError("{} lacks a #CHROM header.".format(path))

    return header_lines, column_header, records


def soft_filter_vcf(
    raw_vcf,
    plain_output,
    record_qc_output,
    expected_sample,
    reference_id,
    reference_length,
):
    header_lines, column_header, records = read_vcf_records(raw_vcf)
    header_sample = column_header.split("\t")[-1]

    if header_sample != expected_sample:
        raise ValueError(
            "VCF sample '{}' does not match expected '{}'."
            .format(header_sample, expected_sample)
        )

    for record in records:
        if record["chrom"] != reference_id:
            raise ValueError(
                "VCF record uses contig '{}', expected '{}'."
                .format(record["chrom"], reference_id)
            )

        if record["pos"] < 1 or record["pos"] > reference_length:
            raise ValueError(
                "VCF position {} outside reference length {}."
                .format(record["pos"], reference_length)
            )

    near_non_snp_positions = set()

    for record in records:
        if record["type"] != "SNP":
            near_non_snp_positions.update(
                circular_positions_around_interval(
                    record["pos"],
                    len(record["ref"]),
                    reference_length,
                    snp_gap,
                )
            )

    filter_counts = Counter()
    type_counts = Counter()
    genotype_counts = Counter()
    qc_rows = []
    pass_haploid_alt_snps = 0

    existing_filter_ids = {
        match.group(1)
        for line in header_lines
        for match in [re.match(r"##FILTER=<ID=([^,>]+)", line)]
        if match
    }

    new_filter_lines = [
        '##FILTER=<ID={},Description="{}">'.format(
            filter_id,
            description.replace('"', "'"),
        )
        for filter_id, description in FILTER_DEFINITIONS.items()
        if filter_id not in existing_filter_ids
    ]

    with plain_output.open("w") as output:
        for line in header_lines:
            output.write(line + "\n")

        output.write(
            "##mitochondrial_consensus_filter="
            "QUAL>={};DP>={};ALT_DP>={};ALT_FRACTION>={};"
            "biallelic_SNP;haploid_ALT_GT;SNP_GAP={}\n"
            .format(
                min_variant_qual,
                min_depth,
                min_alt_depth,
                min_alt_fraction,
                snp_gap,
            )
        )

        for line in new_filter_lines:
            output.write(line + "\n")

        output.write(column_header + "\n")

        for record in records:
            fields = list(record["fields"])
            format_values = record["format"]
            gt = format_values.get("GT", ".")
            dp = parse_int(format_values.get("DP"))
            ad = parse_int_list(format_values.get("AD"))
            adf = parse_int_list(format_values.get("ADF"))
            adr = parse_int_list(format_values.get("ADR"))
            sp = parse_int(format_values.get("SP"))
            scr = parse_int(format_values.get("SCR"))

            try:
                qual = (
                    None
                    if record["qual_text"] in ("", ".")
                    else float(record["qual_text"])
                )
            except ValueError:
                raise ValueError(
                    "Nonnumeric QUAL at {}:{}: {!r}."
                    .format(
                        record["chrom"],
                        record["pos"],
                        record["qual_text"],
                    )
                )

            filters = []

            if record["filter_text"] not in ("", ".", "PASS"):
                filters.extend(record["filter_text"].split(";"))

            if qual is None or qual < min_variant_qual:
                filters.append("LowQual")

            if dp is None or dp < min_depth:
                filters.append("LowDepth")

            if gt == "." or "." in gt:
                filters.append("MissingGT")
                genotype_category = "missing"

            elif "/" in gt or "|" in gt:
                filters.append("NonHaploidGT")
                genotype_category = "non_haploid"

            elif gt == "0":
                filters.append("ReferenceGT")
                genotype_category = "haploid_reference"

            elif gt == "1":
                genotype_category = "haploid_alt"

            else:
                filters.append("NonHaploidGT")
                genotype_category = "unexpected_haploid_allele"

            if record["type"] != "SNP":
                filters.append("NonSNP")

            alt_depth = None
            alt_fraction = None

            if record["type"] == "SNP" and gt == "1":
                if ad is None or len(ad) != 2:
                    filters.append("MissingAlleleDepth")
                else:
                    ref_depth, alt_depth = ad
                    allele_depth_total = ref_depth + alt_depth

                    if alt_depth < min_alt_depth:
                        filters.append("LowAltDepth")

                    if allele_depth_total <= 0:
                        filters.append("MissingAlleleDepth")
                    else:
                        alt_fraction = alt_depth / allele_depth_total

                        if alt_fraction < min_alt_fraction:
                            filters.append("LowAltFraction")

                if record["pos"] in near_non_snp_positions:
                    filters.append("SnpNearIndel")

            # Preserve order while removing duplicate FILTER labels.
            seen_filters = set()
            filters = [
                value
                for value in filters
                if not (value in seen_filters or seen_filters.add(value))
            ]

            final_filter = "PASS" if not filters else ";".join(filters)
            fields[6] = final_filter
            output.write("\t".join(fields) + "\n")

            for value in filters:
                filter_counts[value] += 1

            type_counts[record["type"]] += 1
            genotype_counts[genotype_category] += 1

            candidate = (
                final_filter == "PASS"
                and record["type"] == "SNP"
                and gt == "1"
            )

            if candidate:
                pass_haploid_alt_snps += 1

            qc_rows.append({
                "sample": expected_sample,
                "chrom": record["chrom"],
                "position_1based": record["pos"],
                "REF": record["ref"],
                "ALT": ",".join(record["alts"]),
                "variant_type": record["type"],
                "QUAL": record["qual_text"],
                "FILTER": final_filter,
                "GT": gt,
                "DP": "" if dp is None else dp,
                "AD": "" if ad is None else ",".join(map(str, ad)),
                "ADF": "" if adf is None else ",".join(map(str, adf)),
                "ADR": "" if adr is None else ",".join(map(str, adr)),
                "SP": "" if sp is None else sp,
                "SCR": "" if scr is None else scr,
                "ALT_depth": "" if alt_depth is None else alt_depth,
                "ALT_fraction_AD": (
                    "" if alt_fraction is None else alt_fraction
                ),
                "candidate_for_masked_consensus": (
                    "yes" if candidate else "no"
                ),
            })

    with record_qc_output.open("w", newline="") as handle:
        fieldnames = [
            "sample",
            "chrom",
            "position_1based",
            "REF",
            "ALT",
            "variant_type",
            "QUAL",
            "FILTER",
            "GT",
            "DP",
            "AD",
            "ADF",
            "ADR",
            "SP",
            "SCR",
            "ALT_depth",
            "ALT_fraction_AD",
            "candidate_for_masked_consensus",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    return {
        "total_records": len(records),
        "snp_records": type_counts["SNP"],
        "indel_records": type_counts["INDEL"],
        "complex_records": type_counts["COMPLEX"],
        "pass_haploid_alt_snp_records": pass_haploid_alt_snps,
        "filtered_records": sum(
            1 for row in qc_rows if row["FILTER"] != "PASS"
        ),
        "haploid_reference_gt_records": genotype_counts[
            "haploid_reference"
        ],
        "haploid_alt_gt_records": genotype_counts["haploid_alt"],
        "missing_gt_records": genotype_counts["missing"],
        "non_haploid_or_unexpected_gt_records": (
            genotype_counts["non_haploid"]
            + genotype_counts["unexpected_haploid_allele"]
        ),
        "filter_counts": dict(sorted(filter_counts.items())),
    }


reference_id, reference_length = read_reference(reference_path)
sample_ids = read_sample_ids(sample_list_path)
summary_rows = []
failed_samples = []

for sample in sample_ids:
    bam_path = results_dir / sample / bam_subdir / (sample + bam_suffix)
    final_dir = results_dir / sample / vcf_subdir
    final_vcf = final_dir / (sample + final_vcf_suffix)

    row = {
        "sample": sample,
        "status": "failed",
        "message": "",
        "bam_file": str(bam_path),
        "vcf_file": str(final_vcf),
        "reference_id": reference_id,
        "reference_length_nt": reference_length,
        "bcftools_version": bcftools_version,
        "samtools_version": samtools_version,
        "minimum_mapping_quality": min_map_qual,
        "minimum_base_quality": min_base_qual,
        "minimum_variant_QUAL": min_variant_qual,
        "minimum_depth": min_depth,
        "minimum_alt_depth": min_alt_depth,
        "minimum_alt_fraction": min_alt_fraction,
        "snp_gap_around_non_snp": snp_gap,
        "count_orphans": "yes" if count_orphans else "no",
        "full_BAQ": "yes" if full_baq else "no",
        "observed_max_quality_filtered_depth": "",
        "mpileup_max_depth": "",
        "total_variant_records": "",
        "snp_records": "",
        "indel_records": "",
        "complex_records": "",
        "pass_haploid_alt_snp_records": "",
        "filtered_records": "",
        "haploid_reference_gt_records": "",
        "haploid_alt_gt_records": "",
        "missing_gt_records": "",
        "non_haploid_or_unexpected_gt_records": "",
        "filter_counts": "",
        "record_qc_file": "",
        "bcftools_stats_file": "",
        "commands_log": "",
    }

    if overwrite and final_dir.exists():
        shutil.rmtree(str(final_dir))

    temp_dir = None

    try:
        if not bam_path.is_file():
            raise FileNotFoundError(
                "Missing mapped BAM: {}".format(bam_path)
            )

        quickcheck = subprocess.run(
            ["samtools", "quickcheck", "-v", str(bam_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        if quickcheck.returncode != 0:
            raise ValueError(
                "samtools quickcheck failed: {} {}"
                .format(
                    quickcheck.stdout.strip(),
                    quickcheck.stderr.strip(),
                )
            )

        bam_index = Path(str(bam_path) + ".bai")

        if not bam_index.is_file():
            raise FileNotFoundError(
                "Missing BAM index: {}".format(bam_index)
            )

        header = command_output(
            ["samtools", "view", "-H", str(bam_path)]
        )
        rg_count, bam_samples, markdup_present = parse_bam_header_sample(
            header
        )

        if rg_count < 1:
            raise ValueError(
                "BAM contains no @RG records. Rebuild with map_all_ref.slurm."
            )

        if bam_samples != {sample}:
            raise ValueError(
                "BAM SM values {} do not exactly match sample '{}'."
                .format(sorted(bam_samples), sample)
            )

        if not markdup_present:
            raise ValueError(
                "BAM header does not document duplicate marking. "
                "Rebuild with the revised map_all_ref.slurm."
            )

        idxstats = command_output(
            ["samtools", "idxstats", str(bam_path)]
        )
        reference_seen = False

        for line in idxstats.splitlines():
            if not line.strip():
                continue

            chrom, length_text, mapped_text, unmapped_text = line.split("\t")

            if chrom == "*":
                continue

            length = int(length_text)
            mapped = int(mapped_text)

            if chrom == reference_id:
                reference_seen = True

                if length != reference_length:
                    raise ValueError(
                        "BAM reference length {} differs from FASTA length {}."
                        .format(length, reference_length)
                    )

            elif mapped > 0:
                raise ValueError(
                    "BAM contains mapped reads on unexpected contig '{}'."
                    .format(chrom)
                )

        if not reference_seen:
            raise ValueError(
                "Reference contig '{}' is absent from BAM."
                .format(reference_id)
            )

        depths = observed_depths(bam_path, reference_id)

        if len(depths) != reference_length:
            raise ValueError(
                "samtools depth returned {} positions; expected {}."
                .format(len(depths), reference_length)
            )

        observed_max_depth = max(depths)
        mpileup_max_depth = max(
            mpileup_depth_floor,
            int(
                math.ceil(
                    observed_max_depth * mpileup_depth_multiplier
                )
            ) + mpileup_depth_padding,
        )

        row["observed_max_quality_filtered_depth"] = observed_max_depth
        row["mpileup_max_depth"] = mpileup_max_depth

        sample_root = results_dir / sample
        sample_root.mkdir(parents=True, exist_ok=True)

        temp_dir = Path(tempfile.mkdtemp(
            prefix=".vcf_call.",
            dir=str(sample_root),
        ))

        raw_vcf = temp_dir / (sample + ".raw.normalized.vcf.gz")
        plain_filtered_vcf = temp_dir / (sample + ".softfiltered.vcf")
        temp_final_vcf = temp_dir / final_vcf.name
        record_qc = temp_dir / (sample + ".variant_qc.tsv")
        stats_file = temp_dir / (sample + ".bcftools_stats.txt")
        commands_log = temp_dir / (sample + ".vcf_commands.log")
        sample_ploidy = temp_dir / (sample + ".ploidy.txt")
        sort_temp = temp_dir / "bcftools_sort"
        sort_temp.mkdir()

        sample_ploidy.write_text("{}\t1\n".format(sample))

        mpileup_command = [
            "bcftools", "mpileup",
            "-f", str(reference_path),
            "-r", reference_id,
            "-q", str(min_map_qual),
            "-Q", str(min_base_qual),
            "-d", str(mpileup_max_depth),
            "--threads", str(max(1, threads - 1)),
            "-a",
            (
                "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,"
                "FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,"
                "INFO/SCR"
            ),
        ]

        if count_orphans:
            mpileup_command.append("--count-orphans")

        if full_baq:
            mpileup_command.append("--full-BAQ")

        mpileup_command.extend([
            "-Ou",
            str(bam_path),
        ])

        call_command = [
            "bcftools", "call",
            "--multiallelic-caller",
            "--variants-only",
            "--samples-file", str(sample_ploidy),
            "-Ou",
        ]

        norm_command = [
            "bcftools", "norm",
            "-f", str(reference_path),
            "-c", "e",
            "-m", "-any",
            "-Ou",
        ]

        sort_command = [
            "bcftools", "sort",
            "-T", str(sort_temp),
            "-Oz",
            "-o", str(raw_vcf),
        ]

        with commands_log.open("w") as log:
            run_pipeline(
                [
                    mpileup_command,
                    call_command,
                    norm_command,
                    sort_command,
                ],
                log,
            )

            run_command([
                "bcftools", "index",
                "-f",
                str(raw_vcf),
            ], log)

        raw_samples = [
            line.strip()
            for line in command_output([
                "bcftools", "query",
                "-l",
                str(raw_vcf),
            ]).splitlines()
            if line.strip()
        ]

        if raw_samples != [sample]:
            raise ValueError(
                "Raw VCF samples {} do not equal expected [{}]."
                .format(raw_samples, sample)
            )

        metrics = soft_filter_vcf(
            raw_vcf,
            plain_filtered_vcf,
            record_qc,
            sample,
            reference_id,
            reference_length,
        )

        with commands_log.open("a") as log:
            run_command([
                "bcftools", "view",
                "-Oz",
                "--threads", str(max(1, threads - 1)),
                "-o", str(temp_final_vcf),
                str(plain_filtered_vcf),
            ], log)

            run_command([
                "bcftools", "index",
                "-f",
                str(temp_final_vcf),
            ], log)

            run_command([
                "bcftools", "stats",
                str(temp_final_vcf),
            ], log, stdout_path=stats_file)

        final_samples = [
            line.strip()
            for line in command_output([
                "bcftools", "query",
                "-l",
                str(temp_final_vcf),
            ]).splitlines()
            if line.strip()
        ]

        if final_samples != [sample]:
            raise ValueError(
                "Final VCF samples {} do not equal expected [{}]."
                .format(final_samples, sample)
            )

        # Confirm all final VCF genotypes remain haploid or missing.
        genotype_output = command_output([
            "bcftools", "query",
            "-f", r"[%GT\n]",
            str(temp_final_vcf),
        ])

        invalid_genotypes = sorted({
            gt.strip()
            for gt in genotype_output.splitlines()
            if gt.strip()
            and gt.strip() != "."
            and ("/" in gt or "|" in gt)
        })

        if invalid_genotypes:
            raise ValueError(
                "Diploid/non-haploid genotypes found despite haploid calling: {}."
                .format(",".join(invalid_genotypes))
            )

        # Remove transient files that are not part of the audit output, then
        # replace the entire per-sample VCF directory in one filesystem move.
        if plain_filtered_vcf.exists():
            plain_filtered_vcf.unlink()

        if sort_temp.exists():
            shutil.rmtree(str(sort_temp))

        if final_dir.exists():
            if not overwrite:
                raise FileExistsError(
                    "VCF directory exists and OVERWRITE=false: {}."
                    .format(final_dir)
                )
            shutil.rmtree(str(final_dir))

        shutil.move(str(temp_dir), str(final_dir))
        temp_dir = None

        moved_final = final_dir / temp_final_vcf.name
        moved_index = final_dir / (temp_final_vcf.name + ".csi")
        moved_record_qc = final_dir / record_qc.name
        moved_stats = final_dir / stats_file.name
        moved_commands = final_dir / commands_log.name

        for required in (
            moved_final,
            moved_index,
            moved_record_qc,
            moved_stats,
            moved_commands,
        ):
            if not required.is_file() or required.stat().st_size == 0:
                raise ValueError(
                    "Required finalized VCF output missing or empty: {}."
                    .format(required)
                )

        row.update({
            "status": "passed",
            "message": (
                "haploid normalized audit VCF created; PASS records are "
                "candidate consensus SNPs"
            ),
            **{
                key: value
                for key, value in metrics.items()
                if key != "filter_counts"
            },
            "filter_counts": json.dumps(
                metrics["filter_counts"],
                sort_keys=True,
                separators=(",", ":"),
            ),
            "record_qc_file": str(moved_record_qc),
            "bcftools_stats_file": str(moved_stats),
            "commands_log": str(moved_commands),
        })

    except Exception as exc:
        row["message"] = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    finally:
        if temp_dir is not None:
            shutil.rmtree(str(temp_dir), ignore_errors=True)

    summary_rows.append(row)

summary_path.parent.mkdir(parents=True, exist_ok=True)

with summary_path.open("w", newline="") as handle:
    fieldnames = list(summary_rows[0].keys())
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

print("Reference ID:          {}".format(reference_id))
print("Reference length:      {} bp".format(reference_length))
print("Samples listed:        {}".format(len(sample_ids)))
print("Samples passed:        {}".format(passed_count))
print("Samples failed:        {}".format(len(failed_samples)))
print("Minimum MAPQ:          {}".format(min_map_qual))
print("Minimum base quality:  {}".format(min_base_qual))
print("Minimum variant QUAL:  {}".format(min_variant_qual))
print("Minimum ALT fraction:  {:.3f}".format(min_alt_fraction))
print("Global summary:        {}".format(summary_path))
print("Failed samples:        {}".format(failed_path))

if failed_samples and fail_job_on_error:
    print(
        "ERROR: variant-calling failures: {}."
        .format(", ".join(failed_samples)),
        file=sys.stderr,
    )
    sys.exit(2)
PY
PYTHON_STATUS=$?
set -e

############################################
# FINAL CHECKS
############################################
[[ -s "$SUMMARY" ]] || {
    echo "ERROR: VCF summary missing or empty: $SUMMARY" >&2
    exit 1
}

[[ -f "$FAILED_LIST" ]] || {
    echo "ERROR: failed-sample list was not created: $FAILED_LIST" >&2
    exit 1
}

PASSED=$(awk -F'\t' 'NR>1 && $2=="passed" {n++} END{print n+0}' "$SUMMARY")
FAILED=$(awk -F'\t' 'NR>1 && $2=="failed" {n++} END{print n+0}' "$SUMMARY")

echo
echo "Variant-calling summary:"
echo "  Passed: $PASSED"
echo "  Failed: $FAILED"
echo "  Report: $SUMMARY"

if [[ "$PYTHON_STATUS" -ne 0 ]]; then
    exit "$PYTHON_STATUS"
fi

echo
echo "All samples produced validated haploid mitochondrial audit VCFs."
echo "Next required step: make_masked_consensus.slurm"
