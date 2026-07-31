#!/usr/bin/env bash
#SBATCH --job-name=extract_masked_cds
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/extract_masked_cds_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/extract_masked_cds_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Extract the 13 mitochondrial protein-coding genes from each fixed-coordinate
# masked consensus using the coding-only BED12 produced by cds_bed.slurm.
#
# Output contract for the revised concat_cds.sh:
#   results/SAMPLE/genes_masked/SAMPLE_13PCG_masked.fasta
#
# The output contains exactly 13 records:
#   >SAMPLE|ND1
#   >SAMPLE|ND2
#   ...
#   >SAMPLE|CYTB
#
# This script processes ingroup samples only. Outgroups are handled later by
# the revised outgroup_from_gb.sh using each outgroup's own GenBank annotation.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"

# Outputs created by the revised cds_bed.slurm.
BED="results/cds_coords.bed"
CDS_METADATA="results/cds_coords.tsv"

RESULTS_DIR="results"

# Keep these names synchronized with the revised make_masked_consensus.slurm.
CONSENSUS_SUBDIR="consensus_masked"
CONSENSUS_SUFFIX="_mt_masked.fasta"

GENE_SUBDIR="genes_masked"
GENE_SUFFIX="_13PCG_masked.fasta"
GENE_QC_SUFFIX="_13PCG_qc.tsv"

SUMMARY="$RESULTS_DIR/extract_masked_cds_summary.tsv"
FAILED_LIST="$RESULTS_DIR/extract_masked_cds_failed_samples.txt"

# This step is inexpensive. Rebuild and revalidate outputs on every run rather
# than trusting files left by an older version of the workflow.
OVERWRITE="true"

# Fail the SLURM job after writing all reports if any sample fails validation.
# This prevents downstream scripts from treating a partial extraction as final.
FAIL_JOB_ON_SAMPLE_ERROR="true"

EXPECTED_GENES=(
    "ND1"
    "ND2"
    "COX1"
    "COX2"
    "ATP8"
    "ATP6"
    "COX3"
    "ND3"
    "ND4L"
    "ND4"
    "ND5"
    "ND6"
    "CYTB"
)

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found in Conda environment: $CONDA_ENV" >&2
    exit 1
}

python3 - <<'PY_CHECK'
import sys

if sys.version_info < (3, 6):
    raise SystemExit("ERROR: Python 3.6 or newer is required.")

try:
    import Bio
except ImportError as exc:
    raise SystemExit(
        "ERROR: Biopython is required in the active environment. "
        "Install it with: conda install -c conda-forge biopython"
    ) from exc
PY_CHECK

############################################
# RESOLVE AND CHECK PATHS
############################################
SAMPLE_LIST=$(realpath "$SAMPLE_LIST")
REF=$(realpath "$REF")
BED=$(realpath "$BED")
CDS_METADATA=$(realpath "$CDS_METADATA")
RESULTS_DIR=$(realpath -m "$RESULTS_DIR")
SUMMARY=$(realpath -m "$SUMMARY")
FAILED_LIST=$(realpath -m "$FAILED_LIST")

for file in "$SAMPLE_LIST" "$REF" "$BED" "$CDS_METADATA"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required input is missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$RESULTS_DIR"

EXPECTED_GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

############################################
# EXTRACT AND VALIDATE ALL SAMPLES
############################################
set +e
python3 - \
    "$SAMPLE_LIST" \
    "$REF" \
    "$BED" \
    "$CDS_METADATA" \
    "$RESULTS_DIR" \
    "$CONSENSUS_SUBDIR" \
    "$CONSENSUS_SUFFIX" \
    "$GENE_SUBDIR" \
    "$GENE_SUFFIX" \
    "$GENE_QC_SUFFIX" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$OVERWRITE" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" \
    "$EXPECTED_GENES_TEXT" <<'PY'
import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


sample_list = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
bed_path = Path(sys.argv[3])
metadata_path = Path(sys.argv[4])
results_dir = Path(sys.argv[5])
consensus_subdir = sys.argv[6]
consensus_suffix = sys.argv[7]
gene_subdir = sys.argv[8]
gene_suffix = sys.argv[9]
gene_qc_suffix = sys.argv[10]
summary_path = Path(sys.argv[11])
failed_path = Path(sys.argv[12])
overwrite = sys.argv[13].lower() == "true"
fail_job_on_error = sys.argv[14].lower() == "true"
expected_genes = [x for x in sys.argv[15].split(",") if x]

allowed_bases = set("ACGTRYSWKMBDHVN")
unambiguous_bases = set("ACGT")


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )
    return records[0]


def read_sample_ids(path):
    sample_ids = []
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
            sample_ids.append(sample)

    if not sample_ids:
        raise ValueError("No sample identifiers were found in {}.".format(path))

    return sample_ids


def parse_bed12(path):
    features = []
    seen = set()
    chrom_values = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 12:
                raise ValueError(
                    "{} line {} has {} BED fields; BED12 is required."
                    .format(path, line_number, len(fields))
                )

            chrom = fields[0]
            chrom_start = int(fields[1])
            chrom_end = int(fields[2])
            gene = fields[3].upper()
            strand = fields[5]
            block_count = int(fields[9])
            block_sizes = [
                int(x) for x in fields[10].rstrip(",").split(",") if x
            ]
            block_starts = [
                int(x) for x in fields[11].rstrip(",").split(",") if x
            ]

            if gene in seen:
                raise ValueError(
                    "Duplicate gene '{}' in {}.".format(gene, path)
                )

            if strand not in ("+", "-"):
                raise ValueError(
                    "{} line {} has invalid strand '{}'."
                    .format(path, line_number, strand)
                )

            if chrom_start < 0 or chrom_end <= chrom_start:
                raise ValueError(
                    "{} line {} has invalid interval {}-{}."
                    .format(path, line_number, chrom_start, chrom_end)
                )

            if (
                block_count != len(block_sizes)
                or block_count != len(block_starts)
            ):
                raise ValueError(
                    "{} line {} has inconsistent BED12 block fields."
                    .format(path, line_number)
                )

            blocks = []

            for size, relative_start in zip(block_sizes, block_starts):
                start = chrom_start + relative_start
                end = start + size

                if size <= 0:
                    raise ValueError(
                        "{} line {} contains a non-positive block size."
                        .format(path, line_number)
                    )

                if start < chrom_start or end > chrom_end:
                    raise ValueError(
                        "{} line {} contains a block outside its span."
                        .format(path, line_number)
                    )

                blocks.append((start, end))

            # BED12 blocks must be in ascending genomic order.
            if blocks != sorted(blocks):
                raise ValueError(
                    "{} line {} has unsorted BED12 blocks."
                    .format(path, line_number)
                )

            chrom_values.add(chrom)
            seen.add(gene)
            features.append({
                "gene": gene,
                "chrom": chrom,
                "strand": strand,
                "blocks": blocks,
                "bed_length": sum(end - start for start, end in blocks),
            })

    if len(chrom_values) != 1:
        raise ValueError(
            "{} contains multiple chromosome identifiers: {}."
            .format(path, ", ".join(sorted(chrom_values)))
        )

    observed_genes = [feature["gene"] for feature in features]

    if observed_genes != expected_genes:
        raise ValueError(
            "Gene order/content in {} is incorrect.\nExpected: {}\nObserved: {}"
            .format(
                path,
                ", ".join(expected_genes),
                ", ".join(observed_genes),
            )
        )

    return features, next(iter(chrom_values))


def read_metadata(path):
    metadata = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "gene",
            "strand",
            "coding_nt_retained",
            "transl_table",
            "translation_validated",
        }

        missing_columns = required - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                "{} is missing required columns: {}."
                .format(path, ", ".join(sorted(missing_columns)))
            )

        for row in reader:
            gene = row["gene"].upper()

            if gene in metadata:
                raise ValueError(
                    "Duplicate gene '{}' in {}.".format(gene, path)
                )

            metadata[gene] = {
                "strand": row["strand"],
                "coding_nt_retained": int(row["coding_nt_retained"]),
                "transl_table": int(row["transl_table"]),
                "translation_validated": row["translation_validated"].lower(),
            }

    if list(metadata) != expected_genes:
        raise ValueError(
            "Gene order/content in {} is incorrect.\nExpected: {}\nObserved: {}"
            .format(
                path,
                ", ".join(expected_genes),
                ", ".join(metadata),
            )
        )

    return metadata


def extract_feature(sequence, feature):
    pieces = [
        sequence[start:end]
        for start, end in feature["blocks"]
    ]
    extracted = "".join(pieces)

    if feature["strand"] == "-":
        extracted = str(Seq(extracted).reverse_complement())

    return extracted


def codon_stop_positions(sequence, translation_table):
    stops = []

    for start in range(0, len(sequence), 3):
        codon = sequence[start:start + 3]

        # A definite stop is reportable even when other codons contain N or
        # IUPAC ambiguity symbols.
        if set(codon).issubset(unambiguous_bases):
            amino_acid = str(Seq(codon).translate(table=translation_table))
            if amino_acid == "*":
                stops.append(start // 3 + 1)

    return stops


reference_record = read_single_fasta(reference_path)
reference_sequence = str(reference_record.seq).upper().replace("U", "T")
reference_length = len(reference_sequence)

if not reference_sequence:
    raise ValueError("Reference sequence is empty.")

invalid_reference = sorted(set(reference_sequence) - allowed_bases)
if invalid_reference:
    raise ValueError(
        "Reference contains invalid symbols: {}."
        .format(",".join(invalid_reference))
    )

features, bed_chrom = parse_bed12(bed_path)
metadata = read_metadata(metadata_path)

# The BED chromosome must match the reference FASTA identifier. Sample
# consensus headers may differ because extraction is coordinate-based.
if bed_chrom != reference_record.id:
    raise ValueError(
        "BED chromosome '{}' does not match reference FASTA ID '{}'."
        .format(bed_chrom, reference_record.id)
    )

for feature in features:
    gene = feature["gene"]
    meta = metadata[gene]

    if feature["strand"] != meta["strand"]:
        raise ValueError(
            "Strand disagreement for {} between BED ({}) and metadata ({})."
            .format(gene, feature["strand"], meta["strand"])
        )

    if feature["bed_length"] != meta["coding_nt_retained"]:
        raise ValueError(
            "Length disagreement for {}: BED={} nt, metadata={} nt."
            .format(
                gene,
                feature["bed_length"],
                meta["coding_nt_retained"],
            )
        )

    if feature["bed_length"] % 3 != 0:
        raise ValueError(
            "{} BED length {} is not divisible by 3."
            .format(gene, feature["bed_length"])
        )

    if meta["translation_validated"] != "yes":
        raise ValueError(
            "{} was not marked as translation-validated in {}."
            .format(gene, metadata_path)
        )

    for start, end in feature["blocks"]:
        if start < 0 or end > reference_length:
            raise ValueError(
                "{} block {}-{} exceeds reference length {}."
                .format(gene, start, end, reference_length)
            )

# Independently validate the extraction coordinates on the reference.
for feature in features:
    gene = feature["gene"]
    sequence = extract_feature(reference_sequence, feature)
    meta = metadata[gene]

    if len(sequence) != meta["coding_nt_retained"]:
        raise ValueError(
            "Reference extraction length mismatch for {}.".format(gene)
        )

    reference_stops = codon_stop_positions(
        sequence,
        meta["transl_table"],
    )

    if reference_stops:
        raise ValueError(
            "Reference extraction for {} contains stop codon(s) at amino-acid "
            "positions: {}."
            .format(gene, ",".join(map(str, reference_stops)))
        )

sample_ids = read_sample_ids(sample_list)
summary_rows = []
failed_samples = []

for sample in sample_ids:
    consensus_path = (
        results_dir
        / sample
        / consensus_subdir
        / (sample + consensus_suffix)
    )
    gene_dir = results_dir / sample / gene_subdir
    output_fasta = gene_dir / (sample + gene_suffix)
    output_qc = gene_dir / (sample + gene_qc_suffix)

    status = "failed"
    message = ""
    total_nt = 0
    total_n = 0
    total_other_ambiguity = 0
    sample_qc_rows = []

    # Prevent stale outputs from an earlier run from surviving a failed rerun.
    if overwrite:
        for old_output in (output_fasta, output_qc):
            if old_output.exists():
                old_output.unlink()

    try:
        if not consensus_path.is_file():
            raise FileNotFoundError(
                "Missing masked consensus: {}".format(consensus_path)
            )

        consensus_record = read_single_fasta(consensus_path)
        consensus_sequence = (
            str(consensus_record.seq).upper().replace("U", "T")
        )

        if len(consensus_sequence) != reference_length:
            raise ValueError(
                "Consensus length {} differs from fixed reference length {}. "
                "This suggests an old indel-shifted consensus or an upstream "
                "coordinate error."
                .format(len(consensus_sequence), reference_length)
            )

        invalid = sorted(set(consensus_sequence) - allowed_bases)
        if invalid:
            raise ValueError(
                "Consensus contains invalid or aligned symbols: {}."
                .format(",".join(invalid))
            )

        output_records = []

        for feature in features:
            gene = feature["gene"]
            meta = metadata[gene]
            sequence = extract_feature(consensus_sequence, feature)

            if len(sequence) != meta["coding_nt_retained"]:
                raise ValueError(
                    "{} extraction length {} differs from expected {}."
                    .format(
                        gene,
                        len(sequence),
                        meta["coding_nt_retained"],
                    )
                )

            if len(sequence) % 3 != 0:
                raise ValueError(
                    "{} extraction length {} is not divisible by 3."
                    .format(gene, len(sequence))
                )

            n_count = sequence.count("N")
            other_ambiguity_count = sum(
                1 for base in sequence
                if base not in unambiguous_bases and base != "N"
            )
            definite_stop_positions = codon_stop_positions(
                sequence,
                meta["transl_table"],
            )

            if definite_stop_positions:
                raise ValueError(
                    "{} contains definite internal stop codon(s) at amino-acid "
                    "positions {} under translation table {}."
                    .format(
                        gene,
                        ",".join(map(str, definite_stop_positions)),
                        meta["transl_table"],
                    )
                )

            output_records.append(
                SeqRecord(
                    Seq(sequence),
                    id="{}|{}".format(sample, gene),
                    description="",
                )
            )

            total_nt += len(sequence)
            total_n += n_count
            total_other_ambiguity += other_ambiguity_count

            sample_qc_rows.append({
                "sample": sample,
                "gene": gene,
                "strand": feature["strand"],
                "length_nt": len(sequence),
                "N_count": n_count,
                "percent_N": 100.0 * n_count / len(sequence),
                "other_ambiguity_count": other_ambiguity_count,
                "definite_internal_stop_count": 0,
                "translation_table": meta["transl_table"],
                "source_consensus": str(consensus_path),
            })

        if len(output_records) != len(expected_genes):
            raise ValueError(
                "Extracted {} genes; expected {}."
                .format(len(output_records), len(expected_genes))
            )

        if [record.id.rsplit("|", 1)[1] for record in output_records] != expected_genes:
            raise ValueError("Extracted gene order changed unexpectedly.")

        gene_dir.mkdir(parents=True, exist_ok=True)

        # Write both files atomically. Old validated outputs remain untouched
        # if this sample fails before replacement.
        temp_dir = Path(tempfile.mkdtemp(
            prefix=".extract_cds.",
            dir=str(gene_dir),
        ))

        try:
            temp_fasta = temp_dir / output_fasta.name
            temp_qc = temp_dir / output_qc.name

            written = SeqIO.write(
                output_records,
                str(temp_fasta),
                "fasta",
            )

            if written != len(expected_genes):
                raise ValueError(
                    "Wrote {} FASTA records; expected {}."
                    .format(written, len(expected_genes))
                )

            reread = list(SeqIO.parse(str(temp_fasta), "fasta"))
            expected_ids = [
                "{}|{}".format(sample, gene)
                for gene in expected_genes
            ]

            if [record.id for record in reread] != expected_ids:
                raise ValueError(
                    "FASTA record identifiers/order changed during writing."
                )

            with temp_qc.open("w", newline="") as handle:
                fieldnames = [
                    "sample",
                    "gene",
                    "strand",
                    "length_nt",
                    "N_count",
                    "percent_N",
                    "other_ambiguity_count",
                    "definite_internal_stop_count",
                    "translation_table",
                    "source_consensus",
                ]
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(sample_qc_rows)

            if not overwrite and (output_fasta.exists() or output_qc.exists()):
                raise FileExistsError(
                    "Output exists and OVERWRITE=false for {}."
                    .format(sample)
                )

            os.replace(str(temp_fasta), str(output_fasta))
            os.replace(str(temp_qc), str(output_qc))

        finally:
            shutil.rmtree(str(temp_dir), ignore_errors=True)

        status = "passed"
        message = "13 genes extracted and validated"

    except Exception as exc:
        message = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    summary_rows.append({
        "sample": sample,
        "status": status,
        "message": message,
        "consensus_file": str(consensus_path),
        "output_fasta": str(output_fasta),
        "output_qc": str(output_qc),
        "total_coding_nt": total_nt,
        "total_N": total_n,
        "percent_N_13PCG": (
            100.0 * total_n / total_nt if total_nt else ""
        ),
        "other_ambiguity_count": total_other_ambiguity,
    })

# Global reports are always written, including after sample-level failures.
summary_path.parent.mkdir(parents=True, exist_ok=True)

with summary_path.open("w", newline="") as handle:
    fieldnames = [
        "sample",
        "status",
        "message",
        "consensus_file",
        "output_fasta",
        "output_qc",
        "total_coding_nt",
        "total_N",
        "percent_N_13PCG",
        "other_ambiguity_count",
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

print("Reference ID:              {}".format(reference_record.id))
print("Reference length:          {} bp".format(reference_length))
print("Expected genes:            {}".format(len(expected_genes)))
print("Samples listed:            {}".format(len(sample_ids)))
print("Samples passed:            {}".format(passed_count))
print("Samples failed:            {}".format(len(failed_samples)))
print("Global summary:            {}".format(summary_path))
print("Failed-sample list:        {}".format(failed_path))

if failed_samples and fail_job_on_error:
    print(
        "ERROR: sample-level failures occurred: {}."
        .format(", ".join(failed_samples)),
        file=sys.stderr,
    )
    sys.exit(2)
PY
PYTHON_STATUS=$?
set -e

############################################
# FINAL SHELL-LEVEL CHECKS
############################################
[[ -s "$SUMMARY" ]] || {
    echo "ERROR: global extraction summary is missing or empty: $SUMMARY" >&2
    exit 1
}

# A zero-byte failure list means all samples passed and is valid.
[[ -f "$FAILED_LIST" ]] || {
    echo "ERROR: failed-sample list was not created: $FAILED_LIST" >&2
    exit 1
}

PASSED=$(awk -F'\t' 'NR>1 && $2=="passed" {n++} END{print n+0}' "$SUMMARY")
FAILED=$(awk -F'\t' 'NR>1 && $2=="failed" {n++} END{print n+0}' "$SUMMARY")
TOTAL=$((PASSED + FAILED))

echo
echo "Extraction summary:"
echo "  Total samples: $TOTAL"
echo "  Passed:        $PASSED"
echo "  Failed:        $FAILED"
echo "  Report:        $SUMMARY"

if [[ "$PYTHON_STATUS" -ne 0 ]]; then
    echo "ERROR: CDS extraction completed with sample-level or global failures." >&2
    exit "$PYTHON_STATUS"
fi

echo
echo "All samples passed 13-PCG extraction and validation."
echo "Outputs follow the contract required by the revised concat_cds.sh."
