#!/usr/bin/env bash
#SBATCH --job-name=masked_cons_all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/masked_cons_all_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/masked_cons_all_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Build publication-grade, coordinate-preserving masked mitochondrial
# consensuses from the existing mapped BAM and haploid VCF for each ingroup.
#
# Rules:
# - Apply only high-quality, biallelic SNPs with haploid ALT genotype.
# - Never insert or delete sequence relative to the reference.
# - Mask all called indels, MNPs, multiallelic/uncertain records, low-quality
#   calls, missing/non-haploid genotypes, and positions below minimum depth.
# - Calculate depth after base-quality and mapping-quality filtering.
# - Count overlapping paired-end read segments once.
# - Require the final consensus to remain exactly reference length.
#
# This is the consensus used by:
# - extract_cds.slurm
# - filter_from_con.sh
# - concat_cons.sh

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
VCF_SUFFIX=".vcf.gz"

OUT_SUBDIR="consensus_masked"
OUT_SUFFIX="_mt_masked.fasta"

SUMMARY="$RESULTS_DIR/masked_consensus_summary.tsv"
FAILED_LIST="$RESULTS_DIR/masked_consensus_failed_samples.txt"

# Retain the historical 3x threshold, but calculate and report coverage at
# 1x, 3x, 5x, and 10x. A separate 5x sensitivity reconstruction is advisable
# if focal placements depend on very incomplete specimens.
MIN_DEPTH=3
MIN_BASE_QUAL=20
MIN_MAP_QUAL=20
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
# BUILD MASKED CONSENSUSES
############################################
set +e
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$RESULTS_DIR" \
    "$BAM_SUBDIR" \
    "$BAM_SUFFIX" \
    "$VCF_SUBDIR" \
    "$VCF_SUFFIX" \
    "$OUT_SUBDIR" \
    "$OUT_SUFFIX" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$MIN_DEPTH" \
    "$MIN_BASE_QUAL" \
    "$MIN_MAP_QUAL" \
    "$MIN_VARIANT_QUAL" \
    "$REQUIRE_VCF_SAMPLE_MATCH" \
    "$OVERWRITE" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" <<'PY'
import csv
import gzip
import math
import os
import shutil
import statistics
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
bam_subdir = sys.argv[4]
bam_suffix = sys.argv[5]
vcf_subdir = sys.argv[6]
vcf_suffix = sys.argv[7]
out_subdir = sys.argv[8]
out_suffix = sys.argv[9]
summary_path = Path(sys.argv[10])
failed_path = Path(sys.argv[11])
min_depth = int(sys.argv[12])
min_base_qual = int(sys.argv[13])
min_map_qual = int(sys.argv[14])
min_variant_qual = float(sys.argv[15])
require_vcf_sample_match = sys.argv[16].lower() == "true"
overwrite = sys.argv[17].lower() == "true"
fail_job_on_error = sys.argv[18].lower() == "true"

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


def query_vcf_records(path, sample):
    output = command_output([
        "bcftools", "query",
        "-s", sample,
        "-f", r"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n",
        str(path),
    ])

    records = []
    seen = set()

    for line_number, line in enumerate(output.splitlines(), start=1):
        if not line.strip():
            continue

        fields = line.split("\t")

        if len(fields) != 7:
            raise ValueError(
                "Unexpected bcftools query output at line {} from {}: {} fields."
                .format(line_number, path, len(fields))
            )

        chrom, pos_text, ref, alt, qual, filt, gt = fields
        key = (chrom, int(pos_text), ref, alt)

        if key in seen:
            raise ValueError(
                "Duplicate normalized VCF record: {}:{} {}>{}."
                .format(chrom, pos_text, ref, alt)
            )

        seen.add(key)
        records.append({
            "key": key,
            "chrom": chrom,
            "pos_1based": int(pos_text),
            "ref": ref,
            "alt": alt,
            "qual": qual,
            "filter": filt,
            "gt": gt,
        })

    return records


def merge_intervals(intervals):
    """Merge zero-based, half-open intervals grouped by chromosome."""
    if not intervals:
        return []

    ordered = sorted(intervals, key=lambda x: (x[0], x[1], x[2]))
    merged = []

    for chrom, start, end in ordered:
        if start < 0 or end <= start:
            raise ValueError(
                "Invalid BED interval: {} {} {}."
                .format(chrom, start, end)
            )

        if (
            merged
            and merged[-1][0] == chrom
            and start <= merged[-1][2]
        ):
            merged[-1] = (
                chrom,
                merged[-1][1],
                max(merged[-1][2], end),
            )
        else:
            merged.append((chrom, start, end))

    return merged


def write_bed(path, intervals):
    with path.open("w") as handle:
        for chrom, start, end in intervals:
            handle.write("{}\t{}\t{}\n".format(chrom, start, end))


def interval_length(intervals):
    return sum(end - start for _, start, end in intervals)


reference_record = read_single_fasta(reference_path)
reference_sequence = str(reference_record.seq).upper().replace("U", "T")
reference_id = reference_record.id
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
    bam_path = results_dir / sample / bam_subdir / (sample + bam_suffix)
    vcf_path = results_dir / sample / vcf_subdir / (sample + vcf_suffix)
    final_dir = results_dir / sample / out_subdir
    final_fasta = final_dir / (sample + out_suffix)

    row = {
        "sample": sample,
        "status": "failed",
        "message": "",
        "bam_file": str(bam_path),
        "vcf_file": str(vcf_path),
        "vcf_sample": "",
        "reference_id": reference_id,
        "reference_length_nt": reference_length,
        "mapped_reads": "",
        "total_variant_records": "",
        "hq_haploid_snp_records_applied": "",
        "uncertain_variant_records_masked": "",
        "minimum_depth": min_depth,
        "minimum_base_quality": min_base_qual,
        "minimum_mapping_quality": min_map_qual,
        "minimum_variant_QUAL": min_variant_qual,
        "mean_depth": "",
        "median_depth": "",
        "maximum_depth": "",
        "percent_sites_ge1x": "",
        "percent_sites_ge3x": "",
        "percent_sites_ge5x": "",
        "percent_sites_ge10x": "",
        "low_depth_positions": "",
        "combined_masked_positions": "",
        "length_nt": "",
        "N_count": "",
        "percent_N": "",
        "consensus_file": str(final_fasta),
    }

    if overwrite and final_dir.exists():
        shutil.rmtree(str(final_dir))

    temp_dir = None

    try:
        if not bam_path.is_file():
            raise FileNotFoundError(
                "Missing BAM: {}".format(bam_path)
            )

        if not vcf_path.is_file():
            raise FileNotFoundError(
                "Missing VCF: {}".format(vcf_path)
            )

        quickcheck = subprocess.run(
            ["samtools", "quickcheck", "-v", str(bam_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        if quickcheck.returncode != 0:
            raise ValueError(
                "samtools quickcheck failed for {}: {} {}"
                .format(
                    bam_path,
                    quickcheck.stdout.strip(),
                    quickcheck.stderr.strip(),
                )
            )

        if not Path(str(bam_path) + ".bai").is_file():
            completed = subprocess.run(
                ["samtools", "index", str(bam_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            if completed.returncode != 0:
                raise RuntimeError(
                    "samtools index failed for {}: {}"
                    .format(bam_path, completed.stderr.strip())
                )

        idxstats = command_output(
            ["samtools", "idxstats", str(bam_path)]
        )

        mapped_by_contig = {}

        for line in idxstats.splitlines():
            if not line.strip():
                continue

            chrom, length_text, mapped_text, unmapped_text = line.split("\t")
            if chrom == "*":
                continue

            mapped_by_contig[chrom] = {
                "length": int(length_text),
                "mapped": int(mapped_text),
                "unmapped": int(unmapped_text),
            }

        if reference_id not in mapped_by_contig:
            raise ValueError(
                "Reference contig '{}' is absent from BAM index statistics."
                .format(reference_id)
            )

        if mapped_by_contig[reference_id]["length"] != reference_length:
            raise ValueError(
                "BAM reference length {} differs from FASTA length {}."
                .format(
                    mapped_by_contig[reference_id]["length"],
                    reference_length,
                )
            )

        other_mapped = {
            chrom: values["mapped"]
            for chrom, values in mapped_by_contig.items()
            if chrom != reference_id and values["mapped"] > 0
        }

        if other_mapped:
            raise ValueError(
                "Mapped BAM contains reads on unexpected contigs: {}."
                .format(
                    ", ".join(
                        "{}={}".format(chrom, count)
                        for chrom, count in sorted(other_mapped.items())
                    )
                )
            )

        row["mapped_reads"] = mapped_by_contig[reference_id]["mapped"]

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
        row["vcf_sample"] = vcf_sample

        if require_vcf_sample_match and vcf_sample != sample:
            raise ValueError(
                "VCF sample name '{}' does not match expected sample '{}'."
                .format(vcf_sample, sample)
            )

        sample_root = results_dir / sample
        sample_root.mkdir(parents=True, exist_ok=True)

        temp_dir = Path(tempfile.mkdtemp(
            prefix=".consensus_masked.",
            dir=str(sample_root),
        ))

        variants_dir = temp_dir / "variants_qc"
        variants_dir.mkdir()

        normalized_vcf = variants_dir / (sample + ".normalized.vcf.gz")
        hq_snp_vcf = variants_dir / (sample + ".hq_haploid_snps.vcf.gz")
        all_variant_tsv = variants_dir / (sample + ".all_variants.tsv")
        hq_snp_tsv = variants_dir / (sample + ".hq_haploid_snps.tsv")

        depth_gz = temp_dir / (sample + ".depth_qc.tsv.gz")
        low_depth_bed = temp_dir / (sample + ".low_depth.bed")
        uncertain_bed = temp_dir / (sample + ".uncertain_variants.bed")
        combined_mask_bed = temp_dir / (sample + ".combined_mask.bed")
        raw_fasta = temp_dir / (sample + ".masked.reference_header.fasta")
        temp_final = temp_dir / final_fasta.name
        per_sample_qc = temp_dir / (sample + ".masked_consensus_qc.tsv")
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

        all_records = query_vcf_records(normalized_vcf, vcf_sample)
        hq_records = query_vcf_records(hq_snp_vcf, vcf_sample)
        hq_keys = {record["key"] for record in hq_records}

        if len(hq_keys) != len(hq_records):
            raise ValueError("Duplicate records occurred in HQ SNP VCF.")

        with all_variant_tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "chrom", "position_1based", "REF", "ALT",
                "QUAL", "FILTER", "GT", "classification",
            ])

            for record in all_records:
                classification = (
                    "applied_hq_haploid_snp"
                    if record["key"] in hq_keys
                    else "masked_uncertain_or_non_snp"
                )

                writer.writerow([
                    record["chrom"],
                    record["pos_1based"],
                    record["ref"],
                    record["alt"],
                    record["qual"],
                    record["filter"],
                    record["gt"],
                    classification,
                ])

        with hq_snp_tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "chrom", "position_1based", "REF", "ALT",
                "QUAL", "FILTER", "GT",
            ])

            for record in hq_records:
                writer.writerow([
                    record["chrom"],
                    record["pos_1based"],
                    record["ref"],
                    record["alt"],
                    record["qual"],
                    record["filter"],
                    record["gt"],
                ])

        uncertain_intervals = []

        for record in all_records:
            if record["key"] in hq_keys:
                continue

            if record["chrom"] != reference_id:
                raise ValueError(
                    "VCF record contig '{}' differs from reference ID '{}'."
                    .format(record["chrom"], reference_id)
                )

            start = record["pos_1based"] - 1
            end = start + len(record["ref"])
            start = max(0, start)
            end = min(reference_length, end)

            if end <= start:
                raise ValueError(
                    "Invalid uncertain-variant interval for {}:{} {}>{}."
                    .format(
                        record["chrom"],
                        record["pos_1based"],
                        record["ref"],
                        record["alt"],
                    )
                )

            uncertain_intervals.append((reference_id, start, end))

        uncertain_intervals = merge_intervals(uncertain_intervals)
        write_bed(uncertain_bed, uncertain_intervals)

        # Quality-filtered depth. In samtools depth:
        # -q is minimum BASE quality;
        # -Q is minimum MAPPING quality;
        # -s counts overlapping paired-read segments once.
        depth_output = command_output([
            "samtools", "depth",
            "-aa",
            "-s",
            "-q", str(min_base_qual),
            "-Q", str(min_map_qual),
            "-r", reference_id,
            str(bam_path),
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
                    "samtools depth returned unexpected contig '{}'."
                    .format(chrom)
                )

            position = int(position_text)
            depth = int(depth_text)

            if position < 1 or position > reference_length:
                raise ValueError(
                    "Depth position {} outside reference length {}."
                    .format(position, reference_length)
                )

            if depths[position - 1] is not None:
                raise ValueError(
                    "Duplicate depth position {}."
                    .format(position)
                )

            depths[position - 1] = depth

        missing_depth_positions = [
            index + 1
            for index, value in enumerate(depths)
            if value is None
        ]

        if missing_depth_positions:
            raise ValueError(
                "samtools depth did not report {} reference positions; first: {}."
                .format(
                    len(missing_depth_positions),
                    missing_depth_positions[0],
                )
            )

        with gzip.open(str(depth_gz), "wt") as handle:
            handle.write("chrom\tposition_1based\tdepth\n")

            for position, depth in enumerate(depths, start=1):
                handle.write(
                    "{}\t{}\t{}\n".format(reference_id, position, depth)
                )

        low_depth_intervals = []
        interval_start = None

        for index, depth in enumerate(depths):
            if depth < min_depth and interval_start is None:
                interval_start = index

            if depth >= min_depth and interval_start is not None:
                low_depth_intervals.append(
                    (reference_id, interval_start, index)
                )
                interval_start = None

        if interval_start is not None:
            low_depth_intervals.append(
                (reference_id, interval_start, reference_length)
            )

        write_bed(low_depth_bed, low_depth_intervals)

        combined_intervals = merge_intervals(
            low_depth_intervals + uncertain_intervals
        )
        write_bed(combined_mask_bed, combined_intervals)

        with command_log.open("a") as log:
            consensus_command = [
                "bcftools", "consensus",
                "-f", str(reference_path),
                "-s", vcf_sample,
            ]

            if combined_intervals:
                consensus_command.extend([
                    "-m", str(combined_mask_bed),
                    "--mask-with", "N",
                ])

            consensus_command.extend([
                "-o", str(raw_fasta),
                str(hq_snp_vcf),
            ])

            run(consensus_command, log)

        raw_record = read_single_fasta(raw_fasta)
        final_sequence = str(raw_record.seq).upper().replace("U", "T")

        if len(final_sequence) != reference_length:
            raise ValueError(
                "Final consensus length {} differs from reference length {}."
                .format(len(final_sequence), reference_length)
            )

        invalid_final = sorted(set(final_sequence) - allowed_bases)

        if invalid_final:
            raise ValueError(
                "Final consensus contains invalid symbols: {}."
                .format(",".join(invalid_final))
            )

        SeqIO.write(
            [SeqRecord(Seq(final_sequence), id=sample, description="")],
            str(temp_final),
            "fasta",
        )

        # Index before finalizing to verify samtools accepts the FASTA.
        with command_log.open("a") as log:
            run(["samtools", "faidx", str(temp_final)], log)

        n_count = final_sequence.count("N")
        mean_depth = sum(depths) / float(reference_length)
        median_depth = statistics.median(depths)
        maximum_depth = max(depths)

        def percent_at_least(value):
            return 100.0 * sum(depth >= value for depth in depths) / reference_length

        row.update({
            "status": "passed",
            "message": "coordinate-preserving masked consensus validated",
            "total_variant_records": len(all_records),
            "hq_haploid_snp_records_applied": len(hq_records),
            "uncertain_variant_records_masked": (
                len(all_records) - len(hq_records)
            ),
            "mean_depth": mean_depth,
            "median_depth": median_depth,
            "maximum_depth": maximum_depth,
            "percent_sites_ge1x": percent_at_least(1),
            "percent_sites_ge3x": percent_at_least(3),
            "percent_sites_ge5x": percent_at_least(5),
            "percent_sites_ge10x": percent_at_least(10),
            "low_depth_positions": interval_length(low_depth_intervals),
            "combined_masked_positions": interval_length(combined_intervals),
            "length_nt": reference_length,
            "N_count": n_count,
            "percent_N": 100.0 * n_count / reference_length,
        })

        with per_sample_qc.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(row.keys()),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerow(row)

        # The whole directory is replaced only after every check passes.
        shutil.move(str(temp_dir), str(final_dir))
        temp_dir = None

        # temp_final was moved with the directory and retains its correct name.
        moved_final = final_dir / temp_final.name
        moved_fai = final_dir / (temp_final.name + ".fai")

        if not moved_final.is_file() or not moved_fai.is_file():
            raise ValueError(
                "Final FASTA or FAI missing after directory finalization."
            )

    except Exception as exc:
        row["message"] = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    finally:
        if temp_dir is not None:
            shutil.rmtree(str(temp_dir), ignore_errors=True)

    summary_rows.append(row)

summary_path.parent.mkdir(parents=True, exist_ok=True)

with summary_path.open("w", newline="") as handle:
    fieldnames = [
        "sample",
        "status",
        "message",
        "bam_file",
        "vcf_file",
        "vcf_sample",
        "reference_id",
        "reference_length_nt",
        "mapped_reads",
        "total_variant_records",
        "hq_haploid_snp_records_applied",
        "uncertain_variant_records_masked",
        "minimum_depth",
        "minimum_base_quality",
        "minimum_mapping_quality",
        "minimum_variant_QUAL",
        "mean_depth",
        "median_depth",
        "maximum_depth",
        "percent_sites_ge1x",
        "percent_sites_ge3x",
        "percent_sites_ge5x",
        "percent_sites_ge10x",
        "low_depth_positions",
        "combined_masked_positions",
        "length_nt",
        "N_count",
        "percent_N",
        "consensus_file",
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

print("Reference ID:          {}".format(reference_id))
print("Reference length:      {} bp".format(reference_length))
print("Samples listed:        {}".format(len(sample_ids)))
print("Samples passed:        {}".format(passed_count))
print("Samples failed:        {}".format(len(failed_samples)))
print("Minimum depth:         {}x".format(min_depth))
print("Minimum base quality:  {}".format(min_base_qual))
print("Minimum map quality:   {}".format(min_map_qual))
print("Minimum variant QUAL:  {}".format(min_variant_qual))
print("Global summary:        {}".format(summary_path))
print("Failed samples:        {}".format(failed_path))

if failed_samples and fail_job_on_error:
    print(
        "ERROR: masked consensus failures: {}."
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
    echo "ERROR: masked consensus summary missing or empty: $SUMMARY" >&2
    exit 1
}

[[ -f "$FAILED_LIST" ]] || {
    echo "ERROR: failed-sample list was not created: $FAILED_LIST" >&2
    exit 1
}

PASSED=$(awk -F'\t' 'NR>1 && $2=="passed" {n++} END{print n+0}' "$SUMMARY")
FAILED=$(awk -F'\t' 'NR>1 && $2=="failed" {n++} END{print n+0}' "$SUMMARY")

echo
echo "Masked consensus summary:"
echo "  Passed: $PASSED"
echo "  Failed: $FAILED"
echo "  Report: $SUMMARY"

if [[ "$PYTHON_STATUS" -ne 0 ]]; then
    exit "$PYTHON_STATUS"
fi

echo
echo "All coordinate-preserving masked consensuses passed validation."
echo "These outputs are ready for extract_cds.slurm and filter_from_con.sh."
