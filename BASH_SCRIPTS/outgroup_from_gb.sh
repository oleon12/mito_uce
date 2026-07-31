#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Prepare complete mitochondrial outgroups for both analysis branches:
#
# 1. Whole-mitogenome branch
#    - validate that each FASTA and GenBank file represent the same circular
#      molecule;
#    - orient and rotate each outgroup to the same homologous circular origin
#      as the reference, using the biological start of tRNA-Phe;
#    - validate the 13-PCG order and strand pattern against the reference;
#    - write:
#        results/OUTGROUP/consensus_masked/OUTGROUP_mt_masked.fasta
#
# 2. Thirteen-PCG branch
#    - extract each gene from the outgroup's OWN GenBank annotation;
#    - retain coding sequence only, excluding annotated terminal stop bases;
#    - validate the annotated translation under vertebrate mitochondrial
#      translation table 2;
#    - write exactly 13 records named OUTGROUP|GENE:
#        results/OUTGROUP/genes_masked/OUTGROUP_13PCG_masked.fasta
#
# The "consensus_masked" and "genes_masked" directory names are retained only
# for compatibility with the ingroup pipeline. Public complete outgroup
# genomes are not reference-based masked consensuses.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

REF="references/S_ludovici_QCAZ_18312.fasta"
REFGB="references/S_ludovici_QCAZ_18312.gb"

OUTGROUP_DIR="raw_data/outgroups"

# Recommended manifest format, one outgroup per line:
#   SAMPLE_ID  FASTA_PATH  GENBANK_PATH
#
# Blank lines and lines beginning with # are ignored. Paths may be absolute or
# relative to WORKDIR. If this manifest is absent, exact basename-matched
# FASTA/GenBank pairs are discovered in OUTGROUP_DIR.
OUTGROUP_MANIFEST="CONFS/outgroup_list.txt"

RESULTS_DIR="results"
OUTGROUP_SAMPLE_LIST="$RESULTS_DIR/outgroup_samples.txt"
SUMMARY="$RESULTS_DIR/outgroup_extraction_summary.tsv"
FAILED_LIST="$RESULTS_DIR/outgroup_extraction_failed.txt"

# Complete mammalian mitochondrial genomes are conventionally linearized at
# the beginning of tRNA-Phe. The reference used here follows that convention.
# Outgroups are oriented and rotated so the biological start of this homologous
# feature is at sequence position 1. We do NOT force downstream genes such as
# ND1 to identical absolute coordinates, because genuine indels in upstream
# rRNA/tRNA regions can shift their coordinates among species.
ORIGIN_ANCHOR_PRODUCT="tRNA-Phe"

OVERWRITE="true"
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

for command_name in python3 samtools; do
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

echo "SAMtools version:"
samtools --version | head -n 1

############################################
# RESOLVE AND CHECK PATHS
############################################
REF=$(realpath "$REF")
REFGB=$(realpath "$REFGB")
OUTGROUP_DIR=$(realpath -m "$OUTGROUP_DIR")
OUTGROUP_MANIFEST=$(realpath -m "$OUTGROUP_MANIFEST")
RESULTS_DIR=$(realpath -m "$RESULTS_DIR")
OUTGROUP_SAMPLE_LIST=$(realpath -m "$OUTGROUP_SAMPLE_LIST")
SUMMARY=$(realpath -m "$SUMMARY")
FAILED_LIST=$(realpath -m "$FAILED_LIST")

for file in "$REF" "$REFGB"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required reference input missing or empty: $file" >&2
        exit 1
    }
done

if [[ ! -s "$OUTGROUP_MANIFEST" && ! -d "$OUTGROUP_DIR" ]]; then
    echo "ERROR: neither a nonempty outgroup manifest nor an outgroup directory exists." >&2
    echo "Manifest:  $OUTGROUP_MANIFEST" >&2
    echo "Directory: $OUTGROUP_DIR" >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR"

GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

############################################
# VALIDATE AND PREPARE OUTGROUPS
############################################
set +e
python3 - \
    "$WORKDIR" \
    "$REF" \
    "$REFGB" \
    "$OUTGROUP_DIR" \
    "$OUTGROUP_MANIFEST" \
    "$RESULTS_DIR" \
    "$OUTGROUP_SAMPLE_LIST" \
    "$SUMMARY" \
    "$FAILED_LIST" \
    "$ORIGIN_ANCHOR_PRODUCT" \
    "$OVERWRITE" \
    "$FAIL_JOB_ON_SAMPLE_ERROR" \
    "$GENES_TEXT" <<'PY'
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord


workdir = Path(sys.argv[1])
reference_fasta_path = Path(sys.argv[2])
reference_genbank_path = Path(sys.argv[3])
outgroup_dir = Path(sys.argv[4])
manifest_path = Path(sys.argv[5])
results_dir = Path(sys.argv[6])
outgroup_sample_list_path = Path(sys.argv[7])
summary_path = Path(sys.argv[8])
failed_path = Path(sys.argv[9])
origin_anchor_product = sys.argv[10]
overwrite = sys.argv[11].lower() == "true"
fail_job_on_error = sys.argv[12].lower() == "true"
expected_genes = [x for x in sys.argv[13].split(",") if x]

allowed_nt = set("ACGTRYSWKMBDHVN")
unambiguous_nt = set("ACGT")
mitochondrial_table = CodonTable.unambiguous_dna_by_id[2]
mitochondrial_start_codons = set(mitochondrial_table.start_codons)


GENE_ALIASES = {
    "COI": "COX1",
    "CO1": "COX1",
    "COXI": "COX1",
    "COII": "COX2",
    "CO2": "COX2",
    "COXII": "COX2",
    "COIII": "COX3",
    "CO3": "COX3",
    "COXIII": "COX3",
    "CYB": "CYTB",
    "COB": "CYTB",
    "ATPASE6": "ATP6",
    "ATPASE8": "ATP8",
    "NAD1": "ND1",
    "NAD2": "ND2",
    "NAD3": "ND3",
    "NAD4": "ND4",
    "NAD4L": "ND4L",
    "NAD5": "ND5",
    "NAD6": "ND6",
    "NADHDEHYDROGENASESUBUNIT1": "ND1",
    "NADHDEHYDROGENASESUBUNIT2": "ND2",
    "NADHDEHYDROGENASESUBUNIT3": "ND3",
    "NADHDEHYDROGENASESUBUNIT4": "ND4",
    "NADHDEHYDROGENASESUBUNIT4L": "ND4L",
    "NADHDEHYDROGENASESUBUNIT5": "ND5",
    "NADHDEHYDROGENASESUBUNIT6": "ND6",
    "CYTOCHROMECOXIDASESUBUNITI": "COX1",
    "CYTOCHROMECOXIDASESUBUNIT1": "COX1",
    "CYTOCHROMECOXIDASESUBUNITII": "COX2",
    "CYTOCHROMECOXIDASESUBUNIT2": "COX2",
    "CYTOCHROMECOXIDASESUBUNITIII": "COX3",
    "CYTOCHROMECOXIDASESUBUNIT3": "COX3",
    "CYTOCHROMEB": "CYTB",
    "ATPSYNTHASEFSUBUNIT6": "ATP6",
    "ATPSYNTHASEFSUBUNIT8": "ATP8",
    "ATPSYNTHASEF0SUBUNIT6": "ATP6",
    "ATPSYNTHASEF0SUBUNIT8": "ATP8",
    "ATPSYNTHASESUBUNIT6": "ATP6",
    "ATPSYNTHASESUBUNIT8": "ATP8",
}


def canonical_gene_name(value):
    normalized = re.sub(r"[^A-Za-z0-9]", "", value.upper())

    if normalized in expected_genes:
        return normalized

    return GENE_ALIASES.get(normalized, normalized)


def open_text(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(str(path), "rt")
    return path.open("r")


def parse_records(path, fmt):
    with open_text(path) as handle:
        return list(SeqIO.parse(handle, fmt))


def read_single_record(path, fmt):
    records = parse_records(path, fmt)

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one {} record; found {}."
            .format(path, fmt, len(records))
        )

    return records[0]


def normalize_sequence(sequence):
    return str(sequence).upper().replace("U", "T")


def validate_sequence(sequence, label):
    if not sequence:
        raise ValueError("{} is empty.".format(label))

    invalid = sorted(set(sequence) - allowed_nt)

    if invalid:
        raise ValueError(
            "{} contains invalid nucleotide symbols: {}."
            .format(label, ",".join(invalid))
        )


def sha256_text(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def safe_sample_name(sample):
    if not re.fullmatch(r"[A-Za-z0-9._-]+", sample):
        raise ValueError(
            "Unsafe outgroup identifier '{}'. Use only letters, numbers, "
            "periods, underscores, and hyphens."
            .format(sample)
        )


def resolve_path(value):
    path = Path(value)
    return path if path.is_absolute() else workdir / path


def discover_inputs():
    rows = []
    seen_samples = set()

    if manifest_path.is_file() and manifest_path.stat().st_size > 0:
        mode = "manifest"

        with manifest_path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()

                if not stripped or stripped.startswith("#"):
                    continue

                fields = stripped.split()

                if len(fields) != 3:
                    raise ValueError(
                        "{} line {} has {} fields; expected "
                        "SAMPLE_ID FASTA_PATH GENBANK_PATH."
                        .format(manifest_path, line_number, len(fields))
                    )

                sample, fasta_value, genbank_value = fields
                safe_sample_name(sample)

                if sample in seen_samples:
                    raise ValueError(
                        "Duplicate outgroup '{}' in {}."
                        .format(sample, manifest_path)
                    )

                seen_samples.add(sample)
                rows.append({
                    "sample": sample,
                    "fasta": resolve_path(fasta_value),
                    "genbank": resolve_path(genbank_value),
                    "input_mode": mode,
                })

    else:
        mode = "basename_discovery"
        genbank_extensions = (".gb", ".gbk", ".genbank")
        fasta_extensions = (
            ".fasta", ".fa", ".fna",
            ".fasta.gz", ".fa.gz", ".fna.gz",
        )

        genbank_files = []

        for extension in genbank_extensions:
            genbank_files.extend(sorted(outgroup_dir.glob("*" + extension)))

        if not genbank_files:
            raise ValueError(
                "No GenBank files found in {} and no manifest was supplied."
                .format(outgroup_dir)
            )

        for genbank_path in genbank_files:
            sample = genbank_path.name

            for extension in genbank_extensions:
                if sample.endswith(extension):
                    sample = sample[:-len(extension)]
                    break

            safe_sample_name(sample)

            if sample in seen_samples:
                raise ValueError(
                    "Duplicate discovered outgroup basename '{}'."
                    .format(sample)
                )

            candidates = [
                outgroup_dir / (sample + extension)
                for extension in fasta_extensions
                if (outgroup_dir / (sample + extension)).is_file()
            ]

            if len(candidates) != 1:
                raise ValueError(
                    "Expected exactly one FASTA for '{}'; found {}: {}."
                    .format(
                        sample,
                        len(candidates),
                        ", ".join(map(str, candidates)) or "<none>",
                    )
                )

            seen_samples.add(sample)
            rows.append({
                "sample": sample,
                "fasta": candidates[0],
                "genbank": genbank_path,
                "input_mode": mode,
            })

    if not rows:
        raise ValueError("No outgroup input pairs were identified.")

    for row in rows:
        for label in ("fasta", "genbank"):
            if not row[label].is_file() or row[label].stat().st_size == 0:
                raise ValueError(
                    "Outgroup '{}' {} file is missing or empty: {}."
                    .format(row["sample"], label, row[label])
                )

    return rows


def circular_relation(fasta_sequence, genbank_sequence):
    if len(fasta_sequence) != len(genbank_sequence):
        raise ValueError(
            "FASTA length {} differs from GenBank length {}."
            .format(len(fasta_sequence), len(genbank_sequence))
        )

    length = len(genbank_sequence)

    if fasta_sequence == genbank_sequence:
        return "exact", 0

    doubled = genbank_sequence + genbank_sequence
    position = doubled.find(fasta_sequence)

    if 0 <= position < length:
        return "circular_rotation_of_genbank", position

    reverse = str(Seq(genbank_sequence).reverse_complement())
    doubled_reverse = reverse + reverse
    position = doubled_reverse.find(fasta_sequence)

    if 0 <= position < length:
        return "reverse_complement_circular_rotation_of_genbank", position

    raise ValueError(
        "FASTA and GenBank do not represent the same sequence under exact, "
        "circular-rotation, or reverse-complement circular equivalence."
    )


def exact_location(feature, gene, source):
    parts = list(feature.location.parts)

    if not parts:
        raise ValueError(
            "{} CDS has no location parts in {}."
            .format(gene, source)
        )

    for part in parts:
        if not isinstance(part.start, ExactPosition):
            raise ValueError(
                "{} has a fuzzy CDS start in {}: {}."
                .format(gene, source, part.start)
            )

        if not isinstance(part.end, ExactPosition):
            raise ValueError(
                "{} has a fuzzy CDS end in {}: {}."
                .format(gene, source, part.end)
            )

    if feature.location.strand not in (1, -1):
        raise ValueError(
            "{} has undefined strand in {}."
            .format(gene, source)
        )


def collect_cds_features(record, source):
    features = {}

    for feature in record.features:
        if feature.type != "CDS":
            continue

        candidate_values = []

        for qualifier in ("gene", "locus_tag", "product"):
            candidate_values.extend(feature.qualifiers.get(qualifier, []))

        gene = None

        for value in candidate_values:
            canonical = canonical_gene_name(value)

            if canonical in expected_genes:
                gene = canonical
                break

        if gene is None:
            continue

        exact_location(feature, gene, source)

        if gene in features:
            raise ValueError(
                "Duplicate CDS annotation for {} in {}."
                .format(gene, source)
            )

        features[gene] = feature

    missing = [gene for gene in expected_genes if gene not in features]

    if missing:
        raise ValueError(
            "{} is missing expected CDS annotations: {}."
            .format(source, ", ".join(missing))
        )

    return features


def biological_start(feature):
    strand = feature.location.strand
    parts = list(feature.location.parts)

    if strand == 1:
        return min(int(part.start) for part in parts)

    if strand == -1:
        return max(int(part.end) for part in parts) - 1

    raise ValueError("Feature strand is undefined.")


def find_origin_anchor(record, source, expected_product):
    """
    Find the unique tRNA-Phe feature used as the homologous circular origin.

    Common annotation variants include product=tRNA-Phe and gene=trnF.
    """
    matches = []
    expected_normalized = re.sub(
        r"[^A-Za-z0-9]",
        "",
        expected_product.upper(),
    )

    for feature in record.features:
        if feature.type.lower() != "trna":
            continue

        values = []

        for qualifier in ("product", "gene", "locus_tag", "note"):
            values.extend(feature.qualifiers.get(qualifier, []))

        normalized = [
            re.sub(r"[^A-Za-z0-9]", "", value.upper())
            for value in values
        ]

        is_match = any(
            value == expected_normalized
            or value.startswith(expected_normalized)
            or value in {"TRNF", "TRNAF", "PHENYLALANINETRNA"}
            for value in normalized
        )

        if is_match:
            exact_location(feature, expected_product, source)
            matches.append(feature)

    if len(matches) != 1:
        raise ValueError(
            "{} must contain exactly one '{}' tRNA origin feature; found {}."
            .format(source, expected_product, len(matches))
        )

    return matches[0]


def transformed_gene_geometry(features, genome_length, reverse_complemented, shift):
    starts = {}
    strands = {}

    for gene in expected_genes:
        feature = features[gene]
        start = biological_start(feature)
        strand = feature.location.strand

        if reverse_complemented:
            start = genome_length - 1 - start
            strand = -strand

        start = (start - shift) % genome_length

        starts[gene] = start
        strands[gene] = strand

    order = sorted(expected_genes, key=lambda gene: (starts[gene], gene))
    return starts, strands, order


def translated_sequence_is_valid(coding_sequence, annotated_translation, gene):
    calculated = str(Seq(coding_sequence).translate(table=2))

    if len(calculated) != len(annotated_translation):
        raise ValueError(
            "{} translated length {} differs from annotation length {}."
            .format(gene, len(calculated), len(annotated_translation))
        )

    for index, (observed, expected) in enumerate(
        zip(calculated, annotated_translation)
    ):
        if observed == expected:
            continue

        # GenBank conceptual translations represent an accepted alternative
        # initiation codon as methionine, whereas ordinary codon translation
        # may retain its non-M amino acid unless CDS semantics are requested.
        if (
            index == 0
            and expected == "M"
            and coding_sequence[:3] in mitochondrial_start_codons
        ):
            continue

        raise ValueError(
            "{} translation mismatch at amino-acid position {}: "
            "calculated {}, annotated {}."
            .format(gene, index + 1, observed, expected)
        )

    if "*" in calculated:
        positions = [
            str(index + 1)
            for index, amino_acid in enumerate(calculated)
            if amino_acid == "*"
        ]
        raise ValueError(
            "{} contains internal stop codon(s) at amino-acid positions {}."
            .format(gene, ",".join(positions))
        )


def extract_coding_genes(record, features, source):
    sequences = {}
    rows = []

    for gene in expected_genes:
        feature = features[gene]

        codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
        translation_table = int(
            feature.qualifiers.get("transl_table", ["2"])[0]
        )
        annotated_translation = (
            feature.qualifiers.get("translation", [""])[0]
            .replace(" ", "")
            .replace("\n", "")
            .upper()
        )

        if codon_start != 1:
            raise ValueError(
                "{} in {} has codon_start={}; a complete outgroup CDS with "
                "codon_start=1 is required."
                .format(gene, source, codon_start)
            )

        if translation_table != 2:
            raise ValueError(
                "{} in {} uses translation table {}; expected vertebrate "
                "mitochondrial table 2."
                .format(gene, source, translation_table)
            )

        if not annotated_translation:
            raise ValueError(
                "{} in {} lacks a /translation qualifier."
                .format(gene, source)
            )

        full_sequence = normalize_sequence(feature.extract(record.seq))
        validate_sequence(full_sequence, "{} {} full CDS".format(source, gene))

        coding_length = len(annotated_translation) * 3
        terminal_bases = len(full_sequence) - coding_length

        if terminal_bases not in (0, 1, 2, 3):
            raise ValueError(
                "{} in {} has full CDS length {} but translation implies {} "
                "coding nucleotides; terminal difference {} is invalid."
                .format(
                    gene,
                    source,
                    len(full_sequence),
                    coding_length,
                    terminal_bases,
                )
            )

        coding_sequence = full_sequence[:coding_length]

        if len(coding_sequence) % 3 != 0:
            raise ValueError(
                "{} coding sequence length {} is not divisible by 3."
                .format(gene, len(coding_sequence))
            )

        translated_sequence_is_valid(
            coding_sequence,
            annotated_translation,
            gene,
        )

        sequences[gene] = coding_sequence
        rows.append({
            "gene": gene,
            "product": feature.qualifiers.get("product", [""])[0],
            "strand": "+" if feature.location.strand == 1 else "-",
            "full_annotated_nt": len(full_sequence),
            "coding_nt_retained": len(coding_sequence),
            "terminal_stop_nt_removed": terminal_bases,
            "protein_aa": len(annotated_translation),
            "codon_start": codon_start,
            "translation_table": translation_table,
            "N_count": coding_sequence.count("N"),
            "percent_N": (
                100.0 * coding_sequence.count("N") / len(coding_sequence)
            ),
            "translation_validated": "yes",
        })

    return sequences, rows


def run_samtools_faidx(path):
    completed = subprocess.run(
        ["samtools", "faidx", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "samtools faidx failed for {}: {}."
            .format(path, completed.stderr.strip())
        )


# Validate the reference sequence and annotation.
reference_fasta_record = read_single_record(reference_fasta_path, "fasta")
reference_genbank_record = read_single_record(
    reference_genbank_path,
    "genbank",
)

reference_fasta_sequence = normalize_sequence(reference_fasta_record.seq)
reference_genbank_sequence = normalize_sequence(reference_genbank_record.seq)

validate_sequence(reference_fasta_sequence, "reference FASTA")
validate_sequence(reference_genbank_sequence, "reference GenBank")

if reference_fasta_sequence != reference_genbank_sequence:
    raise ValueError(
        "Reference FASTA and GenBank sequences must be exactly identical "
        "because reference coordinates are used throughout the pipeline."
    )

reference_features = collect_cds_features(
    reference_genbank_record,
    reference_genbank_path,
)

reference_origin_anchor = find_origin_anchor(
    reference_genbank_record,
    reference_genbank_path,
    origin_anchor_product,
)
reference_origin_start = biological_start(reference_origin_anchor)
reference_origin_strand = reference_origin_anchor.location.strand

# The entire ingroup pipeline uses this reference without rotation. Therefore,
# the reference itself must already follow the declared circular convention.
if reference_origin_start != 0 or reference_origin_strand != 1:
    raise ValueError(
        "Reference origin anchor '{}' must begin at position 1 on the plus "
        "strand; observed start={} strand={}. Re-standardize the reference "
        "before running this pipeline."
        .format(
            origin_anchor_product,
            reference_origin_start + 1,
            reference_origin_strand,
        )
    )

reference_starts, reference_strands, reference_gene_order = (
    transformed_gene_geometry(
        reference_features,
        len(reference_genbank_sequence),
        reverse_complemented=False,
        shift=0,
    )
)

input_rows = discover_inputs()
summary_rows = []
failed_samples = []
successful_samples = []

for input_row in input_rows:
    sample = input_row["sample"]
    fasta_path = input_row["fasta"]
    genbank_path = input_row["genbank"]

    final_consensus_dir = results_dir / sample / "consensus_masked"
    final_gene_dir = results_dir / sample / "genes_masked"
    final_qc_dir = results_dir / sample / "outgroup_qc"

    final_consensus = (
        final_consensus_dir / (sample + "_mt_masked.fasta")
    )
    final_genes = (
        final_gene_dir / (sample + "_13PCG_masked.fasta")
    )
    final_gene_qc = (
        final_gene_dir / (sample + "_13PCG_qc.tsv")
    )
    final_standardization_qc = (
        final_qc_dir / (sample + "_circular_standardization.tsv")
    )

    row = {
        "sample": sample,
        "status": "failed",
        "message": "",
        "input_mode": input_row["input_mode"],
        "source_fasta": str(fasta_path),
        "source_genbank": str(genbank_path),
        "genbank_record_id": "",
        "organism": "",
        "genome_length_nt": "",
        "source_fasta_sha256": "",
        "source_genbank_sequence_sha256": "",
        "fasta_genbank_relation": "",
        "fasta_rotation_left_nt_relative_to_genbank": "",
        "origin_anchor_feature": origin_anchor_product,
        "reference_origin_anchor_start_1based": reference_origin_start + 1,
        "reference_origin_anchor_strand": (
            "+" if reference_origin_strand == 1 else "-"
        ),
        "whole_genome_orientation_transform": "",
        "whole_genome_rotation_left_nt": "",
        "standardized_origin_anchor_start_1based": "",
        "pcg_order_matches_reference": "",
        "pcg_strands_match_reference": "",
        "total_13PCG_nt": "",
        "total_13PCG_N": "",
        "whole_genome_N": "",
        "consensus_file": str(final_consensus),
        "genes_file": str(final_genes),
        "gene_qc_file": str(final_gene_qc),
        "standardization_qc_file": str(final_standardization_qc),
    }

    if overwrite:
        for directory in (
            final_consensus_dir,
            final_gene_dir,
            final_qc_dir,
        ):
            if directory.exists():
                shutil.rmtree(str(directory))

    temp_root = None

    try:
        fasta_record = read_single_record(fasta_path, "fasta")
        genbank_record = read_single_record(genbank_path, "genbank")

        fasta_sequence = normalize_sequence(fasta_record.seq)
        genbank_sequence = normalize_sequence(genbank_record.seq)

        validate_sequence(
            fasta_sequence,
            "{} source FASTA".format(sample),
        )
        validate_sequence(
            genbank_sequence,
            "{} source GenBank".format(sample),
        )

        relation, fasta_rotation = circular_relation(
            fasta_sequence,
            genbank_sequence,
        )

        features = collect_cds_features(genbank_record, genbank_path)
        coding_sequences, gene_qc_rows = extract_coding_genes(
            genbank_record,
            features,
            genbank_path,
        )

        genome_length = len(genbank_sequence)
        outgroup_origin_anchor = find_origin_anchor(
            genbank_record,
            genbank_path,
            origin_anchor_product,
        )
        outgroup_origin_start = biological_start(outgroup_origin_anchor)
        outgroup_origin_strand = outgroup_origin_anchor.location.strand

        reverse_complemented = (
            outgroup_origin_strand != reference_origin_strand
        )

        if reverse_complemented:
            oriented_sequence = str(
                Seq(genbank_sequence).reverse_complement()
            )
            oriented_origin_start = (
                genome_length - 1 - outgroup_origin_start
            )
            orientation_transform = "reverse_complemented"
        else:
            oriented_sequence = genbank_sequence
            oriented_origin_start = outgroup_origin_start
            orientation_transform = "as_genbank"

        # Rotate only to the biological start of the homologous origin feature.
        # Do not force ND1 or any downstream gene to an identical coordinate.
        rotation_shift = oriented_origin_start % genome_length

        standardized_sequence = (
            oriented_sequence[rotation_shift:]
            + oriented_sequence[:rotation_shift]
        )

        standardized_origin_start = (
            oriented_origin_start - rotation_shift
        ) % genome_length

        if standardized_origin_start != 0:
            raise ValueError(
                "Circular standardization failed for {}: '{}' maps to "
                "position {}, expected position 1."
                .format(
                    sample,
                    origin_anchor_product,
                    standardized_origin_start + 1,
                )
            )

        outgroup_starts, outgroup_strands, outgroup_gene_order = (
            transformed_gene_geometry(
                features,
                genome_length,
                reverse_complemented=reverse_complemented,
                shift=rotation_shift,
            )
        )

        order_matches = outgroup_gene_order == reference_gene_order
        strands_match = all(
            outgroup_strands[gene] == reference_strands[gene]
            for gene in expected_genes
        )

        if not order_matches:
            raise ValueError(
                "{} 13-PCG order differs from the reference after circular "
                "standardization.\nReference: {}\nOutgroup: {}"
                .format(
                    sample,
                    ",".join(reference_gene_order),
                    ",".join(outgroup_gene_order),
                )
            )

        if not strands_match:
            mismatches = [
                gene
                for gene in expected_genes
                if outgroup_strands[gene] != reference_strands[gene]
            ]
            raise ValueError(
                "{} has PCG strand differences after orientation: {}."
                .format(sample, ",".join(mismatches))
            )

        total_13pcg_nt = sum(
            len(coding_sequences[gene])
            for gene in expected_genes
        )
        total_13pcg_n = sum(
            coding_sequences[gene].count("N")
            for gene in expected_genes
        )

        sample_root = results_dir / sample
        sample_root.mkdir(parents=True, exist_ok=True)
        temp_root = Path(tempfile.mkdtemp(
            prefix=".outgroup_prepare.",
            dir=str(sample_root),
        ))

        temp_consensus_dir = temp_root / "consensus_masked"
        temp_gene_dir = temp_root / "genes_masked"
        temp_qc_dir = temp_root / "outgroup_qc"

        for directory in (
            temp_consensus_dir,
            temp_gene_dir,
            temp_qc_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        temp_consensus = (
            temp_consensus_dir / final_consensus.name
        )
        temp_genes = temp_gene_dir / final_genes.name
        temp_gene_qc = temp_gene_dir / final_gene_qc.name
        temp_standardization_qc = (
            temp_qc_dir / final_standardization_qc.name
        )

        SeqIO.write(
            [
                SeqRecord(
                    Seq(standardized_sequence),
                    id=sample,
                    description="",
                )
            ],
            str(temp_consensus),
            "fasta",
        )
        run_samtools_faidx(temp_consensus)

        gene_records = [
            SeqRecord(
                Seq(coding_sequences[gene]),
                id="{}|{}".format(sample, gene),
                description="",
            )
            for gene in expected_genes
        ]

        written = SeqIO.write(
            gene_records,
            str(temp_genes),
            "fasta",
        )

        if written != len(expected_genes):
            raise ValueError(
                "Wrote {} outgroup genes; expected {}."
                .format(written, len(expected_genes))
            )

        reread_genes = list(SeqIO.parse(str(temp_genes), "fasta"))
        expected_ids = [
            "{}|{}".format(sample, gene)
            for gene in expected_genes
        ]

        if [record.id for record in reread_genes] != expected_ids:
            raise ValueError(
                "Outgroup gene FASTA identifiers/order changed during writing."
            )

        with temp_gene_qc.open("w", newline="") as handle:
            fieldnames = [
                "gene",
                "product",
                "strand",
                "full_annotated_nt",
                "coding_nt_retained",
                "terminal_stop_nt_removed",
                "protein_aa",
                "codon_start",
                "translation_table",
                "N_count",
                "percent_N",
                "translation_validated",
            ]
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(gene_qc_rows)

        with temp_standardization_qc.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "sample",
                "reference_id",
                "reference_length_nt",
                "outgroup_genbank_id",
                "outgroup_length_nt",
                "origin_anchor_feature",
                "reference_origin_anchor_start_1based",
                "outgroup_origin_anchor_start_original_1based",
                "orientation_transform",
                "rotation_left_nt_after_orientation",
                "standardized_origin_anchor_start_1based",
                "reference_pcg_order",
                "standardized_outgroup_pcg_order",
                "pcg_order_matches_reference",
                "pcg_strands_match_reference",
                "whole_genome_source_type",
            ])
            writer.writerow([
                sample,
                reference_genbank_record.id,
                len(reference_genbank_sequence),
                genbank_record.id,
                genome_length,
                origin_anchor_product,
                reference_origin_start + 1,
                outgroup_origin_start + 1,
                orientation_transform,
                rotation_shift,
                standardized_origin_start + 1,
                ",".join(reference_gene_order),
                ",".join(outgroup_gene_order),
                "yes",
                "yes",
                "public_complete_genome_standardized_not_masked",
            ])

        # Verify final temporary products before replacing directories.
        final_check_record = read_single_record(temp_consensus, "fasta")

        if final_check_record.id != sample:
            raise ValueError(
                "Standardized outgroup FASTA ID changed unexpectedly."
            )

        if normalize_sequence(final_check_record.seq) != standardized_sequence:
            raise ValueError(
                "Standardized outgroup FASTA sequence changed during writing."
            )

        for final_directory in (
            final_consensus_dir,
            final_gene_dir,
            final_qc_dir,
        ):
            if final_directory.exists():
                if not overwrite:
                    raise FileExistsError(
                        "Output exists and OVERWRITE=false: {}."
                        .format(final_directory)
                    )
                shutil.rmtree(str(final_directory))

        shutil.move(
            str(temp_consensus_dir),
            str(final_consensus_dir),
        )
        shutil.move(
            str(temp_gene_dir),
            str(final_gene_dir),
        )
        shutil.move(
            str(temp_qc_dir),
            str(final_qc_dir),
        )

        row.update({
            "status": "passed",
            "message": (
                "13 PCGs validated; complete genome orientation and circular "
                "origin standardized to reference"
            ),
            "genbank_record_id": genbank_record.id,
            "organism": genbank_record.annotations.get("organism", ""),
            "genome_length_nt": genome_length,
            "source_fasta_sha256": sha256_text(fasta_sequence),
            "source_genbank_sequence_sha256": sha256_text(
                genbank_sequence
            ),
            "fasta_genbank_relation": relation,
            "fasta_rotation_left_nt_relative_to_genbank": fasta_rotation,
            "whole_genome_orientation_transform": orientation_transform,
            "whole_genome_rotation_left_nt": rotation_shift,
            "standardized_origin_anchor_start_1based": (
                standardized_origin_start + 1
            ),
            "pcg_order_matches_reference": "yes",
            "pcg_strands_match_reference": "yes",
            "total_13PCG_nt": total_13pcg_nt,
            "total_13PCG_N": total_13pcg_n,
            "whole_genome_N": standardized_sequence.count("N"),
        })

        successful_samples.append(sample)

    except Exception as exc:
        row["message"] = str(exc).replace("\t", " ").replace("\n", " ")
        failed_samples.append(sample)

    finally:
        if temp_root is not None:
            shutil.rmtree(str(temp_root), ignore_errors=True)

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

outgroup_sample_list_path.write_text(
    "".join(sample + "\n" for sample in successful_samples)
)

failed_path.write_text(
    "".join(sample + "\n" for sample in failed_samples)
)

print("Reference FASTA ID:       {}".format(reference_fasta_record.id))
print("Reference GenBank ID:     {}".format(reference_genbank_record.id))
print("Reference length:         {} bp".format(
    len(reference_genbank_sequence)
))
print("Circular origin anchor:   {}".format(origin_anchor_product))
print("Reference anchor start:   {}".format(reference_origin_start + 1))
print("Input pairs:              {}".format(len(input_rows)))
print("Outgroups passed:         {}".format(len(successful_samples)))
print("Outgroups failed:         {}".format(len(failed_samples)))
print("Canonical outgroup list:  {}".format(outgroup_sample_list_path))
print("Summary:                  {}".format(summary_path))
print("Failed list:              {}".format(failed_path))

if input_rows and input_rows[0]["input_mode"] == "basename_discovery":
    print(
        "WARNING: inputs were discovered by basename. For publication "
        "reproducibility, create CONFS/outgroup_list.txt.",
        file=sys.stderr,
    )

if failed_samples and fail_job_on_error:
    print(
        "ERROR: outgroup preparation failures: {}."
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
    echo "ERROR: outgroup summary missing or empty: $SUMMARY" >&2
    exit 1
}

[[ -f "$FAILED_LIST" ]] || {
    echo "ERROR: failed-outgroup list was not created: $FAILED_LIST" >&2
    exit 1
}

[[ -f "$OUTGROUP_SAMPLE_LIST" ]] || {
    echo "ERROR: canonical outgroup list was not created: $OUTGROUP_SAMPLE_LIST" >&2
    exit 1
}

PASSED=$(awk -F'\t' 'NR>1 && $2=="passed" {n++} END{print n+0}' "$SUMMARY")
FAILED=$(awk -F'\t' 'NR>1 && $2=="failed" {n++} END{print n+0}' "$SUMMARY")

echo
echo "Outgroup preparation summary:"
echo "  Passed: $PASSED"
echo "  Failed: $FAILED"
echo "  Report: $SUMMARY"

if [[ "$PASSED" -eq 0 ]]; then
    echo "ERROR: no outgroup passed validation." >&2
    exit 1
fi

if [[ "$PYTHON_STATUS" -ne 0 ]]; then
    exit "$PYTHON_STATUS"
fi

echo
echo "All outgroups passed annotation, translation, gene-order, orientation,"
echo "and circular-origin validation."
