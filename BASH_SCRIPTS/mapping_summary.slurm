#!/usr/bin/env bash
#SBATCH --job-name=summ_mapping
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/summ_mapping_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/summ_mapping_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Produce a publication-auditable mapping and depth summary for every sample
# in CONFS/sample_list.txt.
#
# Depth is computed using the same evidence rules as the revised masked
# consensus:
#   minimum base quality 20
#   minimum mapping quality 20
#   overlapping paired-end segments counted once
#   duplicates excluded by samtools depth defaults
#
# The script also reports depth at the first and last mitochondrial reference
# windows as a heuristic check for the artificial linearization boundary of
# a circular genome.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
RESULTS_DIR="results"

SUMMARY="$RESULTS_DIR/mapping_summary.tsv"
FAILED_LIST="$RESULTS_DIR/mapping_summary_failed_samples.txt"

MIN_BASE_QUAL=20
MIN_MAP_QUAL=20
EDGE_WINDOW=200
EDGE_RATIO_WARNING=0.50

FAIL_JOB_ON_SAMPLE_ERROR="true"

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 samtools; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: $command_name was not found in $CONDA_ENV." >&2
        exit 1
    }
done

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
# SUMMARIZE
############################################
set +e
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$RESULTS_DIR" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$MIN_BASE_QUAL" \
    "$MIN_MAP_QUAL" \
    "$EDGE_WINDOW" \
    "$EDGE_RATIO_WARNING" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" <<'PY'
import csv
import math
import statistics
import subprocess
import sys
from pathlib import Path


sample_list_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
results_dir = Path(sys.argv[3])
summary_path = Path(sys.argv[4])
failed_path = Path(sys.argv[5])
min_base_qual = int(sys.argv[6])
min_map_qual = int(sys.argv[7])
edge_window = int(sys.argv[8])
edge_ratio_warning = float(sys.argv[9])
fail_job_on_error = sys.argv[10].lower() == "true"


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


def count_reads(bam, include_flags=None, exclude_flags=None, min_mapq=None):
    command = ["samtools", "view", "-c"]

    if include_flags is not None:
        command.extend(["-f", str(include_flags)])

    if exclude_flags is not None:
        command.extend(["-F", str(exclude_flags)])

    if min_mapq is not None:
        command.extend(["-q", str(min_mapq)])

    command.append(str(bam))
    return int(command_output(command).strip())


def read_sample_ids(path):
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
                raise ValueError(
                    "Duplicate sample '{}' in {}."
                    .format(sample, path)
                )

            seen.add(sample)
            samples.append(sample)

    if not samples:
        raise ValueError("No samples found in {}.".format(path))

    return samples


def read_reference_fai(path):
    fai = Path(str(path) + ".fai")

    with fai.open() as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    if len(rows) != 1:
        raise ValueError(
            "{} must describe exactly one reference sequence; found {}."
            .format(fai, len(rows))
        )

    return rows[0][0], int(rows[0][1])


def read_mapping_qc(path):
    if not path.is_file():
        return {}

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if len(rows) != 1:
        raise ValueError(
            "{} must contain exactly one data row."
            .format(path)
        )

    return rows[0]


sample_ids = read_sample_ids(sample_list_path)
reference_id, reference_length = read_reference_fai(reference_path)

if edge_window <= 0 or edge_window * 2 >= reference_length:
    raise ValueError(
        "EDGE_WINDOW={} is invalid for reference length {}."
        .format(edge_window, reference_length)
    )

summary_rows = []
failed_samples = []

for sample in sample_ids:
    bam = results_dir / sample / "bam" / (sample + ".mapped.bam")
    mapping_qc_path = (
        results_dir / sample / "bam" / (sample + ".mapping_qc.tsv")
    )

    row = {
        "sample": sample,
        "status": "failed",
        "message": "",
        "bam_file": str(bam),
        "reference_id": reference_id,
        "reference_length_nt": reference_length,
        "read_group_SM_present": "",
        "markdup_program_in_header": "",
        "input_read_pairs": "",
        "input_read_segments": "",
        "primary_mapped_read_segments": "",
        "mapping_rate_percent": "",
        "duplicate_read_segments": "",
        "duplicate_percent_of_mapped": "",
        "nonduplicate_mapped_read_segments": "",
        "mapq20_nonduplicate_read_segments": "",
        "proper_pair_nonduplicate_read_segments": "",
        "mate_unmapped_nonduplicate_read_segments": "",
        "minimum_base_quality": min_base_qual,
        "minimum_mapping_quality": min_map_qual,
        "mean_depth_all_sites": "",
        "median_depth_all_sites": "",
        "maximum_depth": "",
        "covered_positions_ge1x": "",
        "percent_sites_ge1x": "",
        "covered_positions_ge3x": "",
        "percent_sites_ge3x": "",
        "covered_positions_ge5x": "",
        "percent_sites_ge5x": "",
        "covered_positions_ge10x": "",
        "percent_sites_ge10x": "",
        "mean_depth_first_edge_window": "",
        "mean_depth_internal_region": "",
        "mean_depth_last_edge_window": "",
        "minimum_edge_to_internal_ratio": "",
        "circular_boundary_review_flag": "",
    }

    try:
        if not bam.is_file():
            raise FileNotFoundError("Missing BAM: {}".format(bam))

        quickcheck = subprocess.run(
            ["samtools", "quickcheck", "-v", str(bam)],
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

        if not Path(str(bam) + ".bai").is_file():
            completed = subprocess.run(
                ["samtools", "index", str(bam)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            if completed.returncode != 0:
                raise RuntimeError(
                    "samtools index failed: {}"
                    .format(completed.stderr.strip())
                )

        header = command_output(["samtools", "view", "-H", str(bam)])
        rg_sm_present = False
        markdup_present = False

        for line in header.splitlines():
            fields = line.split("\t")

            if fields and fields[0] == "@RG":
                if "SM:{}".format(sample) in fields[1:]:
                    rg_sm_present = True

            if fields and fields[0] == "@PG":
                program_fields = {
                    field.split(":", 1)[0]: field.split(":", 1)[1]
                    for field in fields[1:]
                    if ":" in field
                }

                program_name = program_fields.get("PN", "").lower()
                command_line = program_fields.get("CL", "").lower()

                if "markdup" in program_name or "samtools markdup" in command_line:
                    markdup_present = True

        row["read_group_SM_present"] = "yes" if rg_sm_present else "no"
        row["markdup_program_in_header"] = "yes" if markdup_present else "no"

        idxstats = command_output(["samtools", "idxstats", str(bam)])
        contigs = {}

        for line in idxstats.splitlines():
            if not line.strip():
                continue

            chrom, length_text, mapped_text, unmapped_text = line.split("\t")

            if chrom == "*":
                continue

            contigs[chrom] = (
                int(length_text),
                int(mapped_text),
                int(unmapped_text),
            )

        if reference_id not in contigs:
            raise ValueError(
                "Reference contig '{}' absent from BAM."
                .format(reference_id)
            )

        if contigs[reference_id][0] != reference_length:
            raise ValueError(
                "BAM reference length {} differs from FASTA length {}."
                .format(contigs[reference_id][0], reference_length)
            )

        unexpected_mapped = {
            chrom: values[1]
            for chrom, values in contigs.items()
            if chrom != reference_id and values[1] > 0
        }

        if unexpected_mapped:
            raise ValueError(
                "Reads map to unexpected contigs: {}."
                .format(
                    ", ".join(
                        "{}={}".format(chrom, count)
                        for chrom, count in sorted(unexpected_mapped.items())
                    )
                )
            )

        unmapped = count_reads(bam, include_flags=4)
        secondary = count_reads(bam, include_flags=256)
        supplementary = count_reads(bam, include_flags=2048)

        if unmapped or secondary or supplementary:
            raise ValueError(
                "BAM contains excluded classes: unmapped={}, secondary={}, "
                "supplementary={}."
                .format(unmapped, secondary, supplementary)
            )

        mapped = count_reads(bam)
        duplicates = count_reads(bam, include_flags=1024)
        nonduplicates = count_reads(bam, exclude_flags=1024)
        mapq20_nonduplicates = count_reads(
            bam,
            exclude_flags=1024,
            min_mapq=20,
        )
        proper_pair_nonduplicates = count_reads(
            bam,
            include_flags=2,
            exclude_flags=1024,
        )
        mate_unmapped_nonduplicates = count_reads(
            bam,
            include_flags=8,
            exclude_flags=1024,
        )

        if mapped == 0:
            raise ValueError("BAM contains zero primary mapped reads.")

        row["primary_mapped_read_segments"] = mapped
        row["duplicate_read_segments"] = duplicates
        row["duplicate_percent_of_mapped"] = 100.0 * duplicates / mapped
        row["nonduplicate_mapped_read_segments"] = nonduplicates
        row["mapq20_nonduplicate_read_segments"] = mapq20_nonduplicates
        row["proper_pair_nonduplicate_read_segments"] = proper_pair_nonduplicates
        row["mate_unmapped_nonduplicate_read_segments"] = (
            mate_unmapped_nonduplicates
        )

        prior_qc = read_mapping_qc(mapping_qc_path)

        if prior_qc:
            row["input_read_pairs"] = prior_qc.get("input_read_pairs", "")
            row["input_read_segments"] = prior_qc.get(
                "input_read_segments", ""
            )
            row["mapping_rate_percent"] = prior_qc.get(
                "mapping_rate_percent", ""
            )

        depth_output = command_output([
            "samtools", "depth",
            "-aa",
            "-s",
            "-q", str(min_base_qual),
            "-Q", str(min_map_qual),
            "-r", reference_id,
            str(bam),
        ])

        depths = [None] * reference_length

        for line_number, line in enumerate(depth_output.splitlines(), start=1):
            if not line.strip():
                continue

            fields = line.split("\t")

            if len(fields) != 3:
                raise ValueError(
                    "Unexpected samtools depth output at line {}."
                    .format(line_number)
                )

            chrom, position_text, depth_text = fields

            if chrom != reference_id:
                raise ValueError(
                    "Unexpected depth contig '{}'."
                    .format(chrom)
                )

            position = int(position_text)
            depth = int(depth_text)

            if position < 1 or position > reference_length:
                raise ValueError(
                    "Depth position {} outside reference."
                    .format(position)
                )

            if depths[position - 1] is not None:
                raise ValueError(
                    "Duplicate depth position {}."
                    .format(position)
                )

            depths[position - 1] = depth

        missing_positions = [
            index + 1
            for index, value in enumerate(depths)
            if value is None
        ]

        if missing_positions:
            raise ValueError(
                "{} reference positions absent from depth output; first {}."
                .format(len(missing_positions), missing_positions[0])
            )

        def count_ge(threshold):
            return sum(depth >= threshold for depth in depths)

        def percent(count):
            return 100.0 * count / reference_length

        ge1 = count_ge(1)
        ge3 = count_ge(3)
        ge5 = count_ge(5)
        ge10 = count_ge(10)

        first_edge = depths[:edge_window]
        internal = depths[edge_window:-edge_window]
        last_edge = depths[-edge_window:]

        mean_first = sum(first_edge) / len(first_edge)
        mean_internal = sum(internal) / len(internal)
        mean_last = sum(last_edge) / len(last_edge)

        if mean_internal > 0:
            edge_ratio = min(mean_first, mean_last) / mean_internal
        else:
            edge_ratio = float("nan")

        boundary_flag = (
            "review_possible_linearization_edge"
            if math.isfinite(edge_ratio)
            and edge_ratio < edge_ratio_warning
            else "none"
        )

        row.update({
            "status": "passed",
            "message": (
                "mapping and quality-filtered depth validated"
                if rg_sm_present and markdup_present
                else "validated_with_header_warning"
            ),
            "mean_depth_all_sites": sum(depths) / reference_length,
            "median_depth_all_sites": statistics.median(depths),
            "maximum_depth": max(depths),
            "covered_positions_ge1x": ge1,
            "percent_sites_ge1x": percent(ge1),
            "covered_positions_ge3x": ge3,
            "percent_sites_ge3x": percent(ge3),
            "covered_positions_ge5x": ge5,
            "percent_sites_ge5x": percent(ge5),
            "covered_positions_ge10x": ge10,
            "percent_sites_ge10x": percent(ge10),
            "mean_depth_first_edge_window": mean_first,
            "mean_depth_internal_region": mean_internal,
            "mean_depth_last_edge_window": mean_last,
            "minimum_edge_to_internal_ratio": edge_ratio,
            "circular_boundary_review_flag": boundary_flag,
        })

    except Exception as exc:
        row["message"] = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    summary_rows.append(row)

summary_path.parent.mkdir(parents=True, exist_ok=True)

fieldnames = list(summary_rows[0].keys())

with summary_path.open("w", newline="") as handle:
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
rg_warnings = sum(
    row["status"] == "passed"
    and row["read_group_SM_present"] != "yes"
    for row in summary_rows
)
markdup_warnings = sum(
    row["status"] == "passed"
    and row["markdup_program_in_header"] != "yes"
    for row in summary_rows
)
edge_warnings = sum(
    row["status"] == "passed"
    and row["circular_boundary_review_flag"] != "none"
    for row in summary_rows
)

print("Reference ID:                {}".format(reference_id))
print("Reference length:            {} bp".format(reference_length))
print("Samples listed:              {}".format(len(sample_ids)))
print("Samples passed:              {}".format(passed_count))
print("Samples failed:              {}".format(len(failed_samples)))
print("Missing/read-group warnings: {}".format(rg_warnings))
print("Missing markdup warnings:    {}".format(markdup_warnings))
print("Circular-edge warnings:      {}".format(edge_warnings))
print("Summary:                     {}".format(summary_path))
print("Failed samples:              {}".format(failed_path))

if failed_samples and fail_job_on_error:
    print(
        "ERROR: mapping-summary failures: {}."
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
    echo "ERROR: mapping summary missing or empty: $SUMMARY" >&2
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
echo "Mapping summary completed successfully."
echo "Inspect read-group, duplicate-marking, and circular-boundary warning columns."
